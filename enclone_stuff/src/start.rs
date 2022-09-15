// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// See README for documentation.

use crate::analyze_dref::analyze_donor_ref;
use crate::disintegrate::disintegrate_onesies;
use crate::fcell::filter_by_fcell;
use crate::filter_umi::filter_umi;
use crate::flag_defective::flag_defective;
use crate::inconsistent::test_vdj_gex_inconsistent;
use crate::populate_features::populate_features;
use crate::some_filters::some_filters;
use debruijn::dna_string::DnaString;
use enclone::allele::{find_alleles, sub_alts};
use enclone::graph_filter::graph_filter;
use enclone::info::build_info;
use enclone::join::join_exacts;
use enclone::misc1::{cross_filter, lookup_heavy_chain_reuse};
use enclone::misc2::{check_for_barcode_reuse, find_exact_subclonotypes, search_for_shm_indels};
use enclone::misc3::sort_tig_bc;
use enclone_args::read_json::parse_json_annotations_files;
use enclone_core::barcode_fate::BarcodeFate;
use enclone_core::defs::{AlleleData, CloneInfo, TigData};
use enclone_core::enclone_structs::{EncloneExacts, EncloneIntermediates, EncloneSetup};
use enclone_core::hcomp::heavy_complexity;
use enclone_print::define_mat::{define_mat, setup_define_mat};
use enclone_print::loupe::make_donor_refs;
use equiv::EquivRel;
use io_utils::{fwriteln, open_for_read};
use itertools::Itertools;
use qd::dd;
use std::{
    collections::HashMap,
    env,
    fs::File,
    io::{BufRead, BufWriter, Write},
    time::Instant,
};
use string_utils::{add_commas, TextUtils};
use vector_utils::{bin_member, erase_if, next_diff12_3, sort_sync2, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// This is a copy of stirling2_ratio_table from the stirling_numbers crate, that has been modified
// to use higher precision internal math.  This has also been speeded up, and in the process
// made less readable.

use qd::Double;
use rayon::prelude::*;

pub fn stirling2_ratio_table_double(n_max: usize) -> Vec<Vec<Double>> {
    let mut s = Vec::<Vec<Double>>::new();
    let zero = dd![0.0];
    let one = dd![1.0];
    for n in 0..=n_max {
        s.push(vec![zero; n + 1]);
    }
    s[0][0] = one;
    let mut z = Vec::<Double>::new();
    let mut n2n1 = vec![dd![0.0]; n_max + 1];
    for (n, nn) in n2n1.iter_mut().enumerate().skip(2) {
        *nn = Double::from((n - 2) as u32) / Double::from((n - 1) as u32);
    }
    let mut k1k = vec![dd![0.0]; n_max];
    for (k, kk) in k1k.iter_mut().enumerate().skip(1) {
        *kk = Double::from((k - 1) as u32) / Double::from(k as u32);
    }
    let mut njn = Vec::<(usize, Double)>::new();
    for i in 0..n_max + 1 {
        njn.push((i, dd![0.0]));
    }
    njn.par_iter_mut().for_each(|res| {
        let n = res.0;
        if n >= 1 {
            let mut p = one;
            for j in 1..=n {
                p *= Double::from(j as u32) / Double::from(n as u32);
            }
            res.1 = p;
        }
    });

    // This is the slow part of the function.

    for n in 1..=n_max {
        s[n][0] = zero;
        for k in 1..n - 1 {
            z[k - 1] *= k1k[k];
        }
        if n >= 2 {
            z.push(n2n1[n].powi((n - 1) as i32));
        }
        for k in 1..n {
            let x = z[k - 1]; // = ((k-1)/k)^(n-1)
            s[n][k] = s[n - 1][k] + s[n - 1][k - 1] * x;
        }
        s[n][n] = njn[n].1;
    }
    s
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn main_enclone_start(setup: EncloneSetup) -> Result<EncloneIntermediates, String> {
    let tr = Instant::now();
    let ctl = &setup.ctl;
    let gex_info = &setup.gex_info;
    let refdata = &setup.refdata;
    let is_bcr = setup.is_bcr;
    let to_ref_index = &setup.to_ref_index;

    // ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

    // Flag defective reference sequences.

    let mut log = Vec::<u8>::new();
    let mut broken = Vec::<bool>::new();
    flag_defective(ctl, refdata, &mut log, &mut broken);
    ctl.perf_stats(&tr, "flagging defective references");

    // Parse the json annotations file.

    let tparse = Instant::now();
    let mut tig_bc = Vec::<Vec<TigData>>::new();
    let mut vdj_cells = Vec::<Vec<String>>::new();
    let mut gex_cells = Vec::<Vec<String>>::new();
    let mut gex_cells_specified = Vec::<bool>::new();
    let mut fate = vec![HashMap::<String, BarcodeFate>::new(); ctl.origin_info.n()];
    parse_json_annotations_files(
        ctl,
        &mut tig_bc,
        refdata,
        to_ref_index,
        &mut vdj_cells,
        &mut gex_cells,
        &mut gex_cells_specified,
        &mut fate,
    )?;
    ctl.perf_stats(&tparse, "loading from json");

    // Populate features.

    let tpop = Instant::now();
    let mut fr1_starts = Vec::<usize>::new();
    let mut fr2_starts = Vec::<Option<usize>>::new();
    let mut fr3_starts = Vec::<Option<usize>>::new();
    let mut cdr1_starts = Vec::<Option<usize>>::new();
    let mut cdr2_starts = Vec::<Option<usize>>::new();
    populate_features(
        ctl,
        refdata,
        &broken,
        &mut fr1_starts,
        &mut fr2_starts,
        &mut fr3_starts,
        &mut cdr1_starts,
        &mut cdr2_starts,
        &mut log,
    )?;
    if ctl.gen_opt.require_unbroken_ok {
        return Ok(EncloneIntermediates::default());
    }
    for tigi in &mut tig_bc {
        for x in tigi {
            x.fr1_start = fr1_starts[x.v_ref_id];
            x.fr2_start = fr2_starts[x.v_ref_id];
            x.fr3_start = fr3_starts[x.v_ref_id];
            x.cdr1_start = cdr1_starts[x.v_ref_id];
            x.cdr2_start = cdr2_starts[x.v_ref_id];
        }
    }
    ctl.perf_stats(&tpop, "populating features");

    // Test for no data.

    let tproto = Instant::now();
    if ctl.origin_info.n() == 0 {
        return Err("\nNo TCR or BCR data have been specified.\n".to_string());
    }

    // Search for SHM indels.

    search_for_shm_indels(ctl, &tig_bc);
    if ctl.gen_opt.indels {
        return Ok(EncloneIntermediates::default());
    }

    // Record fate of non-cells.

    if ctl.gen_opt.ncell {
        for tigi in &tig_bc {
            let bc = &tigi[0].barcode;
            let li = tigi[0].dataset_index;
            if !bin_member(&vdj_cells[li], bc) {
                fate[li].insert(bc.clone(), BarcodeFate::NotAsmCell);
            }
        }
    }

    // Filter using light --> heavy graph.

    graph_filter(ctl, &mut tig_bc, ctl.gen_opt.graph, &mut fate);

    // Sort tig_bc.

    sort_tig_bc(ctl, &mut tig_bc, refdata);

    // Cross filter.

    cross_filter(ctl, &mut tig_bc, &mut fate);

    // Look for barcode reuse.

    check_for_barcode_reuse(ctl, &tig_bc)?;
    ctl.perf_stats(&tproto, "in proto stuff");

    // Find exact subclonotypes.

    let mut exact_clonotypes = find_exact_subclonotypes(ctl, &tig_bc, refdata, &mut fate);
    if ctl.gen_opt.utr_con || ctl.gen_opt.con_con {
        return Ok(EncloneIntermediates::default());
    }
    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an initial exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }

    // Test for consistency between VDJ cells and GEX cells.

    test_vdj_gex_inconsistent(ctl, &tig_bc, &exact_clonotypes, &vdj_cells, gex_info)?;

    // Filter out some foursie artifacts.

    let t = Instant::now();
    let mut to_delete = vec![false; exact_clonotypes.len()];
    let mut twosies = Vec::<(&[u8], &[u8])>::new();
    for ex in &exact_clonotypes {
        if ex.share.len() == 2 && (ex.share[0].left ^ ex.share[1].left) && ex.ncells() >= 10 {
            twosies.push((ex.share[0].seq.as_ref(), ex.share[1].seq.as_ref()));
        }
    }
    unique_sort(&mut twosies);
    for (ex, d) in exact_clonotypes.iter().zip(to_delete.iter_mut()) {
        if ex.share.len() == 4 {
            for (i1, s1) in ex.share.iter().enumerate() {
                for s2 in &ex.share[i1 + 1..4] {
                    if s1.left ^ s2.left {
                        let p = (s1.seq.as_ref(), s2.seq.as_ref());
                        if bin_member(&twosies, &p) {
                            *d = true;
                            for clone in &ex.clones {
                                fate[clone[0].dataset_index]
                                    .insert(clone[0].barcode.clone(), BarcodeFate::FoursieKill);
                            }
                        }
                    }
                }
            }
        }
    }
    if ctl.clono_filt_opt_def.weak_foursies {
        erase_if(&mut exact_clonotypes, &to_delete);
    }

    // Filter if MAX_HEAVIES = 1 set.

    if ctl.gen_opt.max_heavies == 1 {
        let mut to_delete = vec![false; exact_clonotypes.len()];
        for i in 0..exact_clonotypes.len() {
            let ex = &exact_clonotypes[i];
            let mut heavies = 0;
            for j in 0..ex.share.len() {
                if ex.share[j].left {
                    heavies += 1;
                }
            }
            if heavies > 1 {
                to_delete[i] = true;
            }
        }
        erase_if(&mut exact_clonotypes, &to_delete);
    }
    ctl.perf_stats(&t, "filtering foursies");

    // Build info about clonotypes.  Note that this edits the V reference sequence to perform
    // an indel in some cases.

    let tinfo = Instant::now();
    let mut info: Vec<CloneInfo> = build_info(refdata, ctl, &mut exact_clonotypes, &mut fate);
    ctl.perf_stats(&tinfo, "building info");

    // Derive consensus sequences for alternate alleles of V segments.  Then create donor
    // reference sequences for Loupe.

    let talt = Instant::now();
    // {(donor, ref id, alt seq, support, is_ref)}:
    let mut alt_refs = Vec::<(usize, usize, DnaString, usize, bool)>::new();
    if !ctl.gen_opt.no_alt_alleles {
        alt_refs = find_alleles(refdata, ctl, &exact_clonotypes);
    }
    ctl.perf_stats(&talt, "finding alt alleles");
    if !ctl.gen_opt.dref_file.is_empty() {
        let f = File::create(&ctl.gen_opt.dref_file);
        if f.is_err() {
            eprintln!(
                "\nError trying to write ctl.gen_opt.dref_file = {}.",
                ctl.gen_opt.dref_file
            );
        }
        let mut f = BufWriter::new(f.unwrap());
        let mut count = 0;
        for i in 0..alt_refs.len() {
            let donor = alt_refs[i].0;
            let ref_id = alt_refs[i].1;
            if i > 0 && (donor != alt_refs[i - 1].0 || ref_id != alt_refs[i - 1].1) {
                count = 0;
            }
            let alt_seq = &alt_refs[i].2;
            fwriteln!(
                f,
                ">{}:{}:{}:{} (reference record id : donor name : allele number : gene name)\n{}",
                refdata.id[ref_id],
                ctl.origin_info.donor_id[donor],
                count + 1,
                refdata.name[ref_id],
                alt_seq.to_string()
            );
            count += 1;
        }
    }
    let tdonor = Instant::now();
    let drefs = make_donor_refs(&alt_refs, refdata);
    ctl.perf_stats(&tdonor, "making donor refs");

    // Analyze donor reference.

    analyze_donor_ref(refdata, ctl, &alt_refs);

    // Update reference sequences for V segments by substituting in alt alleles if better.

    sub_alts(refdata, ctl, &alt_refs, &mut info, &mut exact_clonotypes);

    // Compute to_bc, which maps (dataset_index, clonotype_id) to {barcodes}.
    // This is intended as a replacement for some old code below.

    let tbc = Instant::now();
    let mut to_bc = HashMap::<(usize, usize), Vec<String>>::new();
    for (i, ex) in exact_clonotypes.iter().enumerate() {
        for clone in &ex.clones {
            let x = &clone[0];
            to_bc
                .entry((x.dataset_index, i))
                .or_default()
                .push(x.barcode.clone());
        }
    }
    ctl.perf_stats(&tbc, "computing to_bc");

    // Make stirling ratio table.  Not sure that fixing the size of this is safe.

    let tsr = Instant::now();
    let sr = stirling2_ratio_table_double(3000);
    ctl.perf_stats(&tsr, "computing stirling number table");

    // Compute complexity.

    let tcomp = Instant::now();
    if ctl.join_alg_opt.comp_filt < 1_000_000 {
        let jun = heavy_complexity(refdata, &exact_clonotypes, ctl, &drefs);
        for u in 0..exact_clonotypes.len() {
            let ex = &mut exact_clonotypes[u];
            for m in 0..ex.share.len() {
                if ex.share.len() == 2 && ex.share[m].left {
                    ex.share[m].jun = jun[u].clone();
                }
            }
        }
    }
    ctl.perf_stats(&tcomp, "computing complexity");

    // Form equivalence relation on exact subclonotypes.  We also keep the raw joins, consisting
    // of pairs of info indices, that were originally joined.

    let mut join_info = Vec::<(usize, usize, bool, Vec<u8>)>::new();
    let mut raw_joins = Vec::<(i32, i32)>::new();
    let mut eq: EquivRel = join_exacts(
        is_bcr,
        &to_bc,
        refdata,
        ctl,
        &exact_clonotypes,
        &info,
        &mut join_info,
        &mut raw_joins,
        &sr,
        &drefs,
    );

    // If NWEAK_ONESIES is not specified, disintegrate certain onesie clonotypes into single cell
    // clonotypes.  This requires editing of exact_clonotypes, info, eq, join_info and raw_joins.

    let mut disintegrated = Vec::<bool>::new();
    disintegrate_onesies(
        ctl,
        &mut disintegrated,
        &mut eq,
        &mut exact_clonotypes,
        &mut info,
        &mut join_info,
        &mut raw_joins,
    );

    // Update to_bc.

    let txxx = Instant::now();
    let mut to_bc = HashMap::<(usize, usize), Vec<String>>::new();
    for (i, ex) in exact_clonotypes.iter().enumerate() {
        for clone in &ex.clones {
            let x = &clone[0];
            to_bc
                .entry((x.dataset_index, i))
                .or_default()
                .push(x.barcode.clone());
        }
    }

    // Restructure raw joins.

    raw_joins.sort_unstable();
    let raw_joins = {
        let mut raw_joins2 = vec![Vec::<usize>::new(); info.len()];
        for r in raw_joins {
            raw_joins2[r.0 as usize].push(r.1 as usize);
            raw_joins2[r.1 as usize].push(r.0 as usize);
        }
        raw_joins2
    };

    // Lock info.

    let info = &info;

    // Lookup for heavy chain reuse (special purpose experimental option).

    lookup_heavy_chain_reuse(ctl, &exact_clonotypes, info, &eq);
    if ctl.gen_opt.heavy_chain_reuse {
        return Ok(EncloneIntermediates::default());
    }
    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in a pre-filter exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }
    ctl.perf_stats(&txxx, "in some odds and ends");

    // Filter B cells based on UMI counts.

    let tumi = Instant::now();
    let mut orbits = Vec::<Vec<i32>>::new();
    filter_umi(
        &eq,
        &mut orbits,
        ctl,
        &mut exact_clonotypes,
        info,
        &mut fate,
    );
    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an post-umi-filter exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }

    // Remove cells that are not called cells by GEX or feature barcodes.

    let mut orbits = {
        let mut orbits2 = Vec::<Vec<i32>>::new();
        for mut o in orbits {
            let mut to_deletex = vec![false; o.len()];
            for (&x, dx) in o.iter().zip(to_deletex.iter_mut()) {
                let x: &CloneInfo = &info[x as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                let mut to_delete = vec![false; ex.ncells()];
                for (clone, d) in ex.clones.iter().take(ex.ncells()).zip(to_delete.iter_mut()) {
                    let li = clone[0].dataset_index;
                    let bc = &clone[0].barcode;
                    if ctl.gen_opt.cellranger {
                        if gex_cells_specified[li] && !bin_member(&gex_cells[li], bc) {
                            *d = true;
                            fate[li].insert(bc.clone(), BarcodeFate::NotGexCell);
                        }
                    } else if !ctl.origin_info.gex_path[li].is_empty() {
                        let gbc = &gex_info.gex_cell_barcodes[li];
                        if !bin_member(gbc, bc) {
                            fate[li].insert(bc.clone(), BarcodeFate::NotGexCell);
                            if !ctl.clono_filt_opt_def.ngex {
                                *d = true;
                            }
                        }
                    }
                }
                erase_if(&mut ex.clones, &to_delete);
                if ex.ncells() == 0 {
                    *dx = true;
                }
            }
            erase_if(&mut o, &to_deletex);
            if !o.is_empty() {
                orbits2.push(o);
            }
        }
        orbits2
    };

    // Filter using constraints imposed by FCELL.

    filter_by_fcell(ctl, &mut orbits, info, &mut exact_clonotypes, gex_info)?;
    ctl.perf_stats(&tumi, "umi filtering and such");

    // Run some filters.

    some_filters(
        &mut orbits,
        is_bcr,
        &to_bc,
        &sr,
        ctl,
        &exact_clonotypes,
        info,
        &raw_joins,
        &eq,
        &disintegrated,
        &mut fate,
        refdata,
        &drefs,
    );

    // Break up clonotypes containing a large number of chains. These are
    // very likely to be false merges
    let orbits: Vec<Vec<i32>> = orbits
        .into_iter()
        .flat_map(|orbit| {
            let (od, exacts) = setup_define_mat(&orbit, info);
            let mat = define_mat(
                is_bcr,
                &to_bc,
                &sr,
                ctl,
                &exact_clonotypes,
                &exacts,
                &od,
                info,
                &raw_joins,
                refdata,
                &drefs,
            );
            let num_chains = mat.len();
            if num_chains < ctl.join_alg_opt.split_max_chains {
                vec![orbit]
            } else {
                let exacts_of_chains = mat
                    .iter()
                    .enumerate()
                    .flat_map(|(chain_num, chain_in_exact)| {
                        exacts
                            .iter()
                            .zip_eq(chain_in_exact.iter())
                            .filter_map(move |(e, chain)| chain.map(|_| (e, chain_num)))
                    })
                    .into_group_map()
                    .into_iter()
                    .map(|(k, v)| (v, *k))
                    .into_group_map();

                let mut group_of_exacts = HashMap::new();
                let mut group_num = 0;
                for (chains, chain_exacts) in exacts_of_chains {
                    if chains.len() == 1 {
                        for e in chain_exacts {
                            group_of_exacts.insert(e, group_num);
                            group_num += 1;
                        }
                    } else {
                        for e in chain_exacts {
                            group_of_exacts.insert(e, group_num);
                        }
                        group_num += 1;
                    }
                }

                let mut groups = vec![vec![]; group_num];

                for (_, exact_clonotype_id, val) in &od {
                    groups[group_of_exacts[exact_clonotype_id]].push(*val);
                }

                // To split every subclonotype
                // od
                //     .into_iter()
                //     .group_by(|o| o.1)
                //     .into_iter()
                //     .map(|(_, vals)| vals.map(|v| v.2).collect())
                //     .collect();
                groups
            }
        })
        .collect();

    // Pre evaluate (PRE_EVAL).

    if ctl.gen_opt.pre_eval || ctl.join_alg_opt.basic_h.is_some() {
        //
        // Echo command.

        if ctl.gen_opt.echo {
            let args: Vec<String> = env::args().collect();
            println!("{}", args.iter().format(" "));
        }
        if ctl.gen_opt.echoc {
            let args: Vec<String> = env::args().collect();
            println!("# {}", args.iter().format(" "));
        }

        // Gather exact subclonotypes.

        let mut exacts = Vec::<Vec<usize>>::new();
        {
            let mut results = Vec::<(usize, Vec<usize>)>::new();
            for i in 0..orbits.len() {
                results.push((i, Vec::<usize>::new()));
            }
            results.par_iter_mut().for_each(|res| {
                let i = res.0;
                let o = &orbits[i];
                let mut od = Vec::<(Vec<usize>, usize, i32)>::new();
                for id in o.iter() {
                    let x: &CloneInfo = &info[*id as usize];
                    od.push((x.origin.clone(), x.clonotype_id, *id));
                }
                od.sort();
                let mut j = 0;
                while j < od.len() {
                    let k = next_diff12_3(&od, j as i32) as usize;
                    res.1.push(od[j].1);
                    j = k;
                }
            });
            for r in results {
                exacts.push(r.1);
            }
        }

        // Apply POST_FILTER.

        if !ctl.gen_opt.post_filter.is_empty() {
            let mut post_filter = Vec::<(String, String)>::new();
            let f = open_for_read![&ctl.gen_opt.post_filter];
            for (i, line) in f.lines().enumerate() {
                let s = line.unwrap();
                if i == 0 {
                    assert_eq!(s, "dataset,barcode");
                } else {
                    post_filter.push((s.before(",").to_string(), s.after(",").to_string()));
                }
            }
            post_filter.sort();
            for ex in exact_clonotypes.iter_mut() {
                let mut to_delete = vec![false; ex.ncells()];
                for (clone, d) in ex.clones.iter().zip(to_delete.iter_mut()) {
                    let x = &clone[0];
                    if !bin_member(
                        &post_filter,
                        &(
                            ctl.origin_info.dataset_id[x.dataset_index].clone(),
                            x.barcode.clone(),
                        ),
                    ) {
                        *d = true;
                    }
                }
                erase_if(&mut ex.clones, &to_delete);
            }
            let mut to_delete = vec![false; exacts.len()];
            for (ex, d) in exacts.iter().zip(to_delete.iter_mut()) {
                let mut ncells = 0;
                for &e in ex {
                    ncells += exact_clonotypes[e].ncells();
                }
                if ncells == 0 {
                    *d = true;
                }
            }
            erase_if(&mut exacts, &to_delete);
        }

        // Set up metrics.

        let mut merges2 = 0;
        let mut mixes = 0;
        let mut mixed_clonotypes = 0;
        let mut mixed_clonotype_sizes = 0;
        let mut cells1 = 0;
        let mut clonotypes1 = 0;
        let mut cells2 = 0;
        let mut clonotypes2 = 0;
        let mut cells_by_donor = vec![0_usize; ctl.origin_info.donor_list.len()];

        // Reverse sort clonotypes by size.

        let mut n = vec![0; exacts.len()];
        for i in 0..exacts.len() {
            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];
                n[i] += ex.ncells();
            }
        }
        sort_sync2(&mut n, &mut exacts);
        exacts.reverse();
        n.reverse();

        // Process clonotypes.

        for i in 0..exacts.len() {
            let mut mixed = false;
            clonotypes1 += 1;
            cells1 += n[i];
            if n[i] >= 2 {
                clonotypes2 += 1;
            }
            cells2 += n[i];
            if ctl.gen_opt.pre_eval_show && n[i] > 1 {
                println!("\nclonotype");
            }
            let mut cells_by_donor_this = vec![0; ctl.origin_info.donor_list.len()];
            for j in 0..exacts[i].len() {
                let ex = &exact_clonotypes[exacts[i][j]];
                if ctl.gen_opt.pre_eval_show && n[i] > 1 {
                    let donor = ex.clones[0][0].donor_index;
                    if donor.is_some() {
                        print!("{}   ", ctl.origin_info.donor_list[donor.unwrap()]);
                    }
                    for k in 0..ex.share.len() {
                        if ex.share[k].left {
                            print!(
                                "{},{}\t",
                                refdata.name[ex.share[k].v_ref_id],
                                refdata.name[ex.share[k].j_ref_id]
                            );
                        }
                    }
                    for k in 0..ex.share.len() {
                        if !ex.share[k].left {
                            print!(
                                "{},{}\t",
                                refdata.name[ex.share[k].v_ref_id],
                                refdata.name[ex.share[k].j_ref_id]
                            );
                        }
                    }
                    println!();
                }
                for k in 0..ex.clones.len() {
                    let x = &ex.clones[k][0];
                    if x.donor_index.is_some() {
                        cells_by_donor[x.donor_index.unwrap()] += 1;
                        cells_by_donor_this[x.donor_index.unwrap()] += 1;
                    }
                }
            }
            let mut mixes_this = 0;
            if ctl.origin_info.donor_list.len() > 1 && ctl.clono_filt_opt_def.donor {
                for j1 in 0..exacts[i].len() {
                    let ex1 = &exact_clonotypes[exacts[i][j1]];
                    for j2 in j1..exacts[i].len() {
                        let ex2 = &exact_clonotypes[exacts[i][j2]];
                        for k1 in 0..ex1.clones.len() {
                            let x1 = &ex1.clones[k1][0];
                            for k2 in 0..ex2.clones.len() {
                                if (j1, k1) < (j2, k2) {
                                    let x2 = &ex2.clones[k2][0];
                                    if x1.donor_index.is_some()
                                        && x2.donor_index.is_some()
                                        && x1.donor_index.unwrap() != x2.donor_index.unwrap()
                                    {
                                        mixes_this += 1;
                                        mixes += 1;
                                        mixed = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if ctl.gen_opt.pre_eval_show && exacts[i].len() > 1 {
                println!("mixes = {mixes_this}");
            }
            let mut merges2_this = 0;
            for n in cells_by_donor_this {
                if n > 1 {
                    merges2_this += (n * (n - 1)) / 2;
                    merges2 += (n * (n - 1)) / 2;
                }
            }
            if ctl.gen_opt.pre_eval_show && exacts[i].len() > 1 {
                println!("merges = {merges2_this}");
            }
            if mixed {
                mixed_clonotypes += 1;
                mixed_clonotype_sizes += n[i];
            }
        }
        let mut cross = 0;
        let mut intra = 0;
        for i1 in 0..cells_by_donor.len() {
            if cells_by_donor[i1] > 1 {
                intra += cells_by_donor[i1] * (cells_by_donor[i1] - 1) / 2;
            }
            for i2 in i1 + 1..cells_by_donor.len() {
                cross += cells_by_donor[i1] * cells_by_donor[i2];
            }
        }
        println!("\nnumber of intradonor comparisons = {}", add_commas(intra));
        println!(
            "number of intradonor cell-cell merges (quadratic) = {}",
            add_commas(merges2)
        );
        println!("number of cross-donor comparisons = {}", add_commas(cross));
        println!(
            "number of cross-donor comparisons that mix donors = {}",
            add_commas(mixes)
        );
        let rate = (mixes as f64) * 1_000_000_000.0 / (cross as f64);
        println!("rate of cross donor mixing = {:.2} x 10^-9", rate);
        let bogus = (intra as f64) * (mixes as f64) / (cross as f64);
        println!(
            "estimated number of false intradonor merges = {}",
            add_commas(bogus.round() as usize)
        );
        println!("number of mixed clonotypes = {mixed_clonotypes}");
        println!(
            "percent of non-single-cell mixed clonotypes = {:.2}",
            100.0 * mixed_clonotypes as f64 / clonotypes2 as f64
        );
        println!("sum of mixed clonotype sizes = {mixed_clonotype_sizes}");
        println!("total number of cells in clonotypes = {cells1}");
        println!(
            "mean clonotype size = {:.3}",
            cells1 as f64 / clonotypes1 as f64
        );
        println!(
            "mean non-single-cell clonotype size = {:.3}\n",
            cells2 as f64 / clonotypes2 as f64
        );
        std::process::exit(0);
    }

    // Mark VDJ noncells.

    let tmark = Instant::now();
    if ctl.clono_filt_opt_def.non_cell_mark {
        for ex in exact_clonotypes.iter_mut() {
            for clone in ex.clones.iter_mut() {
                let di = clone[0].dataset_index;
                if !bin_member(&vdj_cells[di], &clone[0].barcode) {
                    clone[0].marked = true;
                }
            }
        }
    }
    ctl.perf_stats(&tmark, "marking vdj noncells");
    if !ctl.gen_opt.trace_barcode.is_empty() {
        for ex in &exact_clonotypes {
            for clone in &ex.clones {
                if clone[0].barcode == ctl.gen_opt.trace_barcode {
                    println!(
                        "\nfound {} in an intermediate exact subclonotype having {} cells",
                        ctl.gen_opt.trace_barcode,
                        ex.ncells(),
                    );
                }
            }
        }
    }
    Ok(EncloneIntermediates {
        setup,
        ex: EncloneExacts {
            to_bc,
            exact_clonotypes,
            raw_joins,
            info: info.to_vec(),
            orbits,
            vdj_cells,
            join_info,
            drefs,
            sr,
            fate,
            is_bcr,
            allele_data: AlleleData {
                alt_refs,
                var_pos: Vec::new(),
                var_bases: Vec::new(),
            },
        },
    })
}
