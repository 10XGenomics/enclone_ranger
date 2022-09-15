// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Delete exact subclonotypes that appear to represent doublets.
//
// THIS FILTER DOESN'T PROPERLY TRACK FATE.

use enclone_core::{
    barcode_fate::BarcodeFate,
    defs::{CloneInfo, EncloneControl, ExactClonotype},
};
use enclone_print::define_mat::{define_mat, setup_define_mat};
use enclone_proto::types::DonorReferenceItem;
use itertools::Itertools;
use qd::Double;
use rayon::prelude::*;
use std::collections::HashMap;
use std::time::Instant;
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, erase_if, next_diff, next_diff1_2, sort_sync2};

pub fn delete_doublets(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    fate: &mut [HashMap<String, BarcodeFate>],
) {
    if ctl.clono_filt_opt_def.doublet {
        let t = Instant::now();
        // Define pure subclonotypes.  To do this we break each clonotype up by chain signature.
        // Note duplication of code with print_clonotypes.rs.  And this is doing some
        // superfluous compute.

        let mut results = Vec::<(usize, Vec<Vec<usize>>)>::new();
        for i in 0..orbits.len() {
            results.push((i, Vec::new()));
        }
        let mut pures = Vec::<Vec<usize>>::new();
        ctl.perf_stats(&t, "doublet filtering setup");
        let t = Instant::now();
        results.par_iter_mut().for_each(|res| {
            let i = res.0;
            let o = orbits[i].clone();
            let (od, mut exacts) = setup_define_mat(&o, info);
            let mat = define_mat(
                is_bcr,
                to_bc,
                sr,
                ctl,
                exact_clonotypes,
                &exacts,
                &od,
                info,
                raw_joins,
                refdata,
                dref,
            );
            let nexacts = mat[0].len();
            let mut priority = Vec::<Vec<bool>>::new();
            for u in 0..nexacts {
                let mut typex = vec![false; mat.len()];
                for col in 0..mat.len() {
                    if mat[col][u].is_some() {
                        typex[col] = true;
                    }
                }
                priority.push(typex.clone());
            }
            sort_sync2(&mut priority, &mut exacts);
            let mut j = 0;
            while j < priority.len() {
                let k = next_diff(&priority, j);
                let mut p = Vec::<usize>::new();
                for &e in &exacts[j..k] {
                    p.push(e);
                }
                res.1.push(p);
                j = k;
            }
        });
        for mut r in results {
            pures.append(&mut r.1);
        }

        // Define the number of cells in each pure subclonotype.

        let mut npure = vec![0; pures.len()];
        for j in 0..pures.len() {
            for id in pures[j].iter() {
                npure[j] += exact_clonotypes[*id].ncells();
            }
        }

        // Find the pairs of pure subclonotypes that share identical CDR3 sequences.

        ctl.perf_stats(&t, "doublet filtering main");
        let t = Instant::now();
        let mut shares = Vec::<(usize, usize)>::new();
        {
            let mut content = Vec::<(&str, usize)>::new();
            for (j, pure) in pures.iter().enumerate() {
                for &id in pure {
                    let ex = &exact_clonotypes[id];
                    for s in &ex.share {
                        content.push((s.cdr3_dna.as_str(), j));
                    }
                }
            }
            content.par_sort();
            content.dedup();
            ctl.perf_stats(&t, "doublet filtering shares setup");
            let t = Instant::now();
            let mut j = 0;
            while j < content.len() {
                let k = next_diff1_2(&content, j as i32) as usize;
                for l1 in j..k {
                    for l2 in l1 + 1..k {
                        shares.push((content[l1].1, content[l2].1));
                        shares.push((content[l2].1, content[l1].1));
                    }
                }
                j = k;
            }
            shares.par_sort();
            shares.dedup();
            ctl.perf_stats(&t, "doublet filtering shares");
        }

        // Find triples of pure subclonotypes in which the first two have no share, but both
        // of the first two share with the third.

        let t = Instant::now();
        const MIN_MULT_DOUBLET: usize = 5;
        let mut trips = Vec::<(usize, usize, usize)>::new();
        {
            let mut us = Vec::<usize>::new();
            let mut vs = Vec::<Vec<usize>>::new();
            let mut j = 0;
            while j < shares.len() {
                // not using next_diff1_2 here because of i32 overflow issue
                let mut k = j + 1;
                loop {
                    if k == shares.len() || shares[k].0 != shares[j].0 {
                        break;
                    }
                    k += 1;
                }
                let u = shares[j].0;
                us.push(u);
                let mut x = Vec::<usize>::new();
                for v in &shares[j..k] {
                    let v = v.1;
                    if MIN_MULT_DOUBLET * npure[u] <= npure[v] {
                        x.push(v);
                    }
                }
                vs.push(x);
                j = k;
            }
            let mut results = Vec::<(usize, Vec<(usize, usize, usize)>)>::new();
            for i in 0..us.len() {
                results.push((i, Vec::new()));
            }
            results.par_iter_mut().for_each(|res| {
                let i = res.0;
                let u = us[i];
                let vs = &vs[i];
                for l1 in 0..vs.len() {
                    for l2 in l1 + 1..vs.len() {
                        let v1 = vs[l1];
                        let v2 = vs[l2];
                        if !bin_member(&shares, &(v1, v2)) {
                            res.1.push((v1, v2, u));
                        }
                    }
                }
            });
            for mut r in results {
                trips.append(&mut r.1);
            }
        }
        ctl.perf_stats(&t, "doublet filtering trips");

        // Delete some of the third members of the triples.

        let t = Instant::now();
        let mut to_delete = vec![false; exact_clonotypes.len()];
        for (v1, v2, v0) in trips {
            let verbose = false;
            if verbose {
                println!("\n{}, {}, {}", v0, v1, v2);
                println!("DELETING");
                for (u, m) in pures[v0].iter().enumerate() {
                    let ex = &exact_clonotypes[*m];
                    let mut cdrs = Vec::<String>::new();
                    for k in 0..ex.share.len() {
                        cdrs.push(ex.share[k].cdr3_aa.clone());
                    }
                    println!("[{}] {}", u + 1, cdrs.iter().format(","));
                }
                println!("USING");
                for (u, m) in pures[v1].iter().enumerate() {
                    let ex = &exact_clonotypes[*m];
                    let mut cdrs = Vec::<String>::new();
                    for k in 0..ex.share.len() {
                        cdrs.push(ex.share[k].cdr3_aa.clone());
                    }
                    println!("[{}] {}", u + 1, cdrs.iter().format(","));
                }
                println!("AND");
                for (u, m) in pures[v2].iter().enumerate() {
                    let ex = &exact_clonotypes[*m];
                    let mut cdrs = Vec::<String>::new();
                    for k in 0..ex.share.len() {
                        cdrs.push(ex.share[k].cdr3_aa.clone());
                    }
                    println!("[{}] {}", u + 1, cdrs.iter().format(","));
                }
            }
            for m in pures[v0].iter() {
                to_delete[*m] = true;
            }
        }
        let mut orbits2 = Vec::<Vec<i32>>::new();
        for o in orbits.iter() {
            let mut o = o.clone();
            let mut del2 = vec![false; o.len()];
            for j in 0..o.len() {
                let id = info[o[j] as usize].clonotype_index;
                if to_delete[id] {
                    del2[j] = true;
                    let x: &CloneInfo = &info[o[j] as usize];
                    let ex = &exact_clonotypes[x.clonotype_index];
                    for k in 0..ex.ncells() {
                        let li = ex.clones[k][0].dataset_index;
                        let bc = &ex.clones[k][0].barcode;
                        fate[li].insert(bc.clone(), BarcodeFate::Doublet);
                    }
                }
            }
            erase_if(&mut o, &del2);
            orbits2.push(o);
        }
        *orbits = orbits2;
        ctl.perf_stats(&t, "doublet filtering tail");
    }
}
