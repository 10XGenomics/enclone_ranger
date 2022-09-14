// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Based on the number of cells in each column, decide which exact subclonotypes
// look like junk.  Preliminary heuristic.

use enclone_core::{
    barcode_fate::BarcodeFate,
    defs::{CloneInfo, EncloneControl, ExactClonotype},
};
use enclone_print::define_mat::{define_mat, setup_define_mat};
use enclone_proto::types::DonorReferenceItem;
use qd::Double;
use rayon::prelude::*;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::{erase_if, next_diff12_3};

pub fn weak_chains(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    fate: &mut [HashMap<String, BarcodeFate>],
    refdata: &RefData,
    dref: &[DonorReferenceItem],
) {
    // Note mat calculation duplicated with print_clonotypes and also doublet detection.

    let mut results = Vec::<(usize, Vec<(usize, String, BarcodeFate)>, Vec<usize>)>::new();
    for i in 0..orbits.len() {
        results.push((i, Vec::new(), Vec::new()));
    }
    results.par_iter_mut().for_each(|res| {
        let i = res.0;
        let o = orbits[i].clone();
        let (od, exacts) = setup_define_mat(&o, info);
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
        let cols = mat.len();
        if cols > 2 {
            let nexacts = exacts.len();
            let mut ncells = vec![0; cols];
            let mut col_entries = vec![Vec::<usize>::new(); cols];
            for (u, &clonotype_id) in exacts.iter().enumerate().take(nexacts) {
                for ((nc, ce), m) in ncells
                    .iter_mut()
                    .zip(col_entries.iter_mut())
                    .zip(mat.iter().take(cols))
                {
                    let mid = m[u];
                    if mid.is_some() {
                        ce.push(u);
                        *nc += exact_clonotypes[clonotype_id].clones.len();
                    }
                }
            }
            let mut total_cells = 0;
            for j in 0..exacts.len() {
                total_cells += exact_clonotypes[exacts[j]].ncells();
            }
            for j in 0..cols {
                if ncells[j] <= 20 && 8 * ncells[j] < total_cells {
                    for d in col_entries[j].iter() {
                        if ctl.clono_filt_opt_def.weak_chains {
                            res.2.push(exacts[*d]);
                        }
                        let ex = &exact_clonotypes[exacts[*d]];
                        for i in 0..ex.ncells() {
                            res.1.push((
                                ex.clones[i][0].dataset_index,
                                ex.clones[i][0].barcode.clone(),
                                BarcodeFate::WeakChains,
                            ));
                        }
                    }
                }
            }
        }
    });
    let mut to_delete = vec![false; exact_clonotypes.len()];
    let mut dels = Vec::<i32>::new();
    for i in 0..results.len() {
        for j in 0..results[i].1.len() {
            fate[results[i].1[j].0].insert(results[i].1[j].1.clone(), results[i].1[j].2.clone());
        }
        for x in results[i].2.iter() {
            to_delete[*x] = true;
        }
    }
    dels.sort_unstable();
    for o in orbits.iter_mut() {
        let del = o
            .iter()
            .map(|&oj| {
                let id = info[oj as usize].clonotype_index;
                to_delete[id]
            })
            .collect::<Vec<_>>();
        erase_if(o, &del);
    }
}
