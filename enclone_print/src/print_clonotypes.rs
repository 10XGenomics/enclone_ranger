// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

use crate::define_mat::define_mat;
use crate::filter::survives_filter;
use crate::finish_table::finish_table;
use crate::gene_scan::gene_scan_test;
use crate::loupe::{loupe_out, make_loupe_clonotype};
use crate::print_utils1::{compute_field_types, extra_args, start_gen};
use crate::print_utils2::row_fill;
use crate::print_utils3::{
    consensus_codon_cdr3, define_column_info, get_extra_parseables, process_complete,
};
use crate::print_utils4::{build_show_aa, compute_bu, compute_some_stats};
use crate::print_utils5::{delete_weaks, vars_and_shares};
use enclone_args::proc_args_check::involves_gex_fb;
use enclone_core::allowed_vars::{CVARS_ALLOWED, CVARS_ALLOWED_PCELL, LVARS_ALLOWED};
use enclone_core::defs::{AlleleData, CloneInfo, ColInfo, EncloneControl, ExactClonotype, GexInfo};
use enclone_core::mammalian_fixed_len::mammalian_fixed_len_peer_groups;
use enclone_core::set_speakers::set_speakers;
use enclone_proto::types::{Clonotype, DonorReferenceItem};
use equiv::EquivRel;
use hdf5::Reader;
use qd::Double;
use rayon::prelude::*;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufWriter;
use string_utils::TextUtils;
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, bin_position, erase_if, next_diff12_3, unique_sort};

// Print clonotypes.  A key challenge here is to define the columns that represent shared
// chains.  This is given below by the code that forms an equivalence relation on the CDR3_AAs.
//
// This code carries out a second function, which is to filter out exact subclonotypes in orbits
// that appear to be junk.  Exactly how these should be reflected in files is TBD.
//
// Some inputs for this section:
// refdata                = reference sequence data
// ctl                    = control parameters
// exact_clonotypes       = the vector of all exact subclonotypes
// info                   = vector of clonotype info
// eq                     = equivalence relation on info

pub fn print_clonotypes(
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    orbits: &[Vec<i32>],
    raw_joins: &[Vec<usize>],
    gex_info: &GexInfo,
    vdj_cells: &[Vec<String>],
    d_readers: &[Option<Reader>],
    ind_readers: &[Option<Reader>],
    h5_data: &[(usize, Vec<u32>, Vec<u32>)],
    pics: &mut Vec<String>,
    exacts: &mut Vec<Vec<usize>>,
    in_center: &mut Vec<bool>,
    rsi: &mut Vec<ColInfo>,
    out_datas: &mut Vec<Vec<HashMap<String, String>>>,
    tests: &mut Vec<usize>,
    controls: &mut Vec<usize>,
    fate: &mut [HashMap<String, &str>],
    allele_data: &AlleleData,
) -> Result<(), String> {
    let lvars = &ctl.clono_print_opt.lvars;

    // Compute extra args.

    let extra_args = extra_args(ctl);

    // Determine if any lvars need gex info.

    let need_gex = {
        lvars.iter().map(String::as_str).any(involves_gex_fb)
            || {
                if ctl.parseable_opt.pout.is_empty() {
                    false
                } else if ctl.parseable_opt.pcols.is_empty() {
                    LVARS_ALLOWED.into_iter().any(involves_gex_fb)
                } else {
                    ctl.parseable_opt
                        .pcols
                        .iter()
                        .map(String::as_str)
                        .any(involves_gex_fb)
                }
            }
            || extra_args.iter().map(String::as_str).any(involves_gex_fb)
    };

    // Define parseable output columns.  The entire machinery for parseable output is controlled
    // by macros that begin with "speak".

    let mut max_chains = 4;
    // This seems like a bug, since rsi is uninitialized upon entry to print_clonotypes.
    for r in rsi.iter() {
        max_chains = max(max_chains, r.mat.len());
    }
    let mut parseable_fields = Vec::<String>::new();
    set_speakers(ctl, &mut parseable_fields, max_chains);
    let pcols_sort = &ctl.parseable_opt.pcols_sort;

    // Identify certain extra parseable variables.  These arise from parameterizable cvars.

    let mut extra_parseables = get_extra_parseables(ctl, pcols_sort);

    // Compute all_vars.

    let rsi_vars = &ctl.clono_print_opt.cvars;
    let mut all_vars = rsi_vars.iter().map(String::as_str).collect::<Vec<_>>();
    for var in CVARS_ALLOWED {
        if !rsi_vars.contains(&var.to_string()) {
            all_vars.push(var);
        }
    }
    for var in CVARS_ALLOWED_PCELL {
        if !rsi_vars.contains(&var.to_string()) {
            all_vars.push(var);
        }
    }
    all_vars.append(&mut extra_parseables);
    for x in extra_args.iter() {
        if !rsi_vars.contains(x) {
            all_vars.push(x.as_str());
        }
    }

    // Test for presence of GEX/FB data.

    let mut have_gex = false;
    for i in 0..ctl.origin_info.gex_path.len() {
        if !ctl.origin_info.gex_path[i].is_empty() {
            have_gex = true;
        }
    }

    // Gather alt_bcs_fields.

    let mut alt_bcs = Vec::<String>::new();
    for li in 0..ctl.origin_info.alt_bc_fields.len() {
        for i in 0..ctl.origin_info.alt_bc_fields[li].len() {
            alt_bcs.push(ctl.origin_info.alt_bc_fields[li][i].0.clone());
        }
    }
    unique_sort(&mut alt_bcs);

    // Compute number of vdj cells that are gex.

    let mut n_vdj_gex = Vec::<usize>::new();
    for (gex, vdj) in gex_info
        .pca
        .iter()
        .zip(vdj_cells.iter())
        .take(ctl.origin_info.n())
    {
        let mut n = 0;
        for y in gex.iter() {
            if bin_member(vdj, y.0) {
                n += 1;
            }
        }
        n_vdj_gex.push(n);
    }

    // Compute peer groups.

    let peer_groups = mammalian_fixed_len_peer_groups(refdata);

    // Traverse the orbits.

    // 0: index in reps
    // 1: vector of clonotype pictures
    // 2: vector of some clonotype info
    //    [parallel to 1]
    // next to last three entries = whitelist contam, denominator for that, low gex count
    // added out_datas (used to be next to last three, now one more)
    let mut results = Vec::<(
        usize,
        Vec<String>,
        Vec<(Vec<usize>, ColInfo)>,
        usize,
        usize,
        usize,
        Vec<Clonotype>,
        Vec<Vec<HashMap<String, String>>>,
        isize,
        Vec<bool>,
        Vec<bool>,
        Vec<(usize, String, &'static str)>,
        Vec<bool>,
        String,
    )>::new();
    for i in 0..orbits.len() {
        results.push((
            i,
            Vec::<String>::new(),
            Vec::<(Vec<usize>, ColInfo)>::new(),
            0,
            0,
            0,
            Vec::<Clonotype>::new(),
            Vec::<Vec<HashMap<String, String>>>::new(),
            0,
            Vec::<bool>::new(),
            Vec::<bool>::new(),
            Vec::new(),
            Vec::new(),
            String::new(),
        ));
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

        // Reconstruct the participating clones.  This is needed because most exact subclonotypes
        // having more than two chains have been split up.
        //
        // Capture these data into parallel data structures, one per exact subclonotype:
        // exacts: the exact subclonotype ids
        // mults:  number of cells [redundant, might remove]

        let mut exacts = Vec::<usize>::new();
        let mut mults = Vec::<usize>::new();
        let mut j = 0;
        let loupe_clonotypes = &mut res.6;
        while j < od.len() {
            let k = next_diff12_3(&od, j as i32) as usize;
            let mut mult = 0_usize;
            for l in j..k {
                let x: &CloneInfo = &info[od[l].2 as usize];
                let m = x.clonotype_index;
                mult = exact_clonotypes[m].clones.len();
            }
            mults.push(mult);
            exacts.push(od[j].1);
            j = k;
        }

        // There are two passes.  On the first pass we only identify the exact subclonotypes that
        // are junk.  On the second pass we remove those and then print the orbit.

        let mut bads = vec![false; exacts.len()];
        let mut stats_pass1 = Vec::<Vec<(String, Vec<String>)>>::new();
        for pass in 1..=2 {
            // Delete weak exact subclonotypes.

            if pass == 2 && !ctl.clono_filt_opt.protect_bads {
                erase_if(&mut mults, &bads);
                erase_if(&mut exacts, &bads);
            }

            // Sort exact subclonotypes.

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
            let priority = exacts
                .iter()
                .enumerate()
                .map(|(u, &exact)| {
                    let typex = mat.iter().map(|col| col[u].is_some()).collect::<Vec<_>>();
                    let clonotype_id = exact;
                    let ex = &exact_clonotypes[clonotype_id];
                    let mut utot0 = 0;
                    if let Some(mid) = mat[0][u] {
                        let ex = &exact_clonotypes[clonotype_id];
                        for j in 0..ex.clones.len() {
                            utot0 += ex.clones[j][mid].umi_count;
                        }
                    }
                    (typex, ex.ncells(), utot0)
                })
                .collect::<Vec<_>>();
            let permutation = permutation::sort(&priority[..]);
            exacts = permutation.apply_slice(&exacts[..]);
            mults = permutation.apply_slice(&mults[..]);
            exacts.reverse();
            mults.reverse();

            // Define a matrix mat[col][ex] which is the column of the exact subclonotype
            // corresponding to the given column col of the clonotype, which may or may not be
            // defined.  Then define other information associated to each chain.  These are
            // reference sequence identifiers, CDR3 start positions, and the like.

            let nexacts = exacts.len();
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
            let mut rsi = define_column_info(ctl, &exacts, exact_clonotypes, &mat, refdata);
            rsi.mat = mat;
            let mat = &rsi.mat;

            // Let n be the total number of cells in this pass.

            let n: usize = mults.iter().sum();

            // Filter.

            let mut in_center = true;
            if pass == 2
                && !survives_filter(
                    &exacts,
                    &rsi,
                    ctl,
                    exact_clonotypes,
                    refdata,
                    gex_info,
                    dref,
                )
            {
                if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                    in_center = false;
                } else {
                    continue;
                }
            }

            // Generate Loupe data.

            if (!ctl.gen_opt.binary.is_empty() || !ctl.gen_opt.proto.is_empty()) && pass == 2 {
                loupe_clonotypes.push(make_loupe_clonotype(
                    exact_clonotypes,
                    &exacts,
                    &rsi,
                    refdata,
                    dref,
                    ctl,
                ));
            }

            // Set up for parseable output.

            let mut out_data = Vec::<HashMap<String, String>>::new();

            // Print the orbit.
            // ◼ An assumption of this code is that a productive pair does not have two contigs
            // ◼ having identical CDR3_AA sequences.  At present this is not enforced by the
            // ◼ assembly stage, so the assumption is violated.  To work around this there are
            // ◼ some unsavory workarounds below.

            let mut mlog = Vec::<u8>::new();
            if n >= ctl.clono_filt_opt.ncells_low
                || ctl.clono_group_opt.asymmetric_center == "from_filters"
            {
                // Start to generate parseable output.

                if pass == 2 {
                    start_gen(
                        ctl,
                        &exacts,
                        exact_clonotypes,
                        &mut out_data,
                        &mut mlog,
                        &extra_args,
                    );
                }

                // Find variant positions.  And some other things.

                let mut vars = Vec::<Vec<usize>>::new();
                let mut vars_amino = Vec::<Vec<usize>>::new();
                let mut shares_amino = Vec::<Vec<usize>>::new();
                let mut ref_diff_pos = Vec::<Vec<Vec<usize>>>::new();
                vars_and_shares(
                    pass,
                    ctl,
                    &exacts,
                    exact_clonotypes,
                    &rsi,
                    refdata,
                    dref,
                    &mut vars,
                    &mut vars_amino,
                    &mut shares_amino,
                    &mut ref_diff_pos,
                    &mut out_data,
                );

                // Mark some weak exact subclonotypes for deletion.

                if pass == 1 {
                    delete_weaks(ctl, &exacts, exact_clonotypes, mat, refdata, &mut bads);
                }

                // Done unless on second pass.  Unless there are bounds or COMPLETE specified
                // or VAR_DEF specified.

                if pass == 1
                    && ctl.clono_filt_opt.bounds.is_empty()
                    && !ctl.gen_opt.complete
                    && ctl.gen_opt.var_def.is_empty()
                {
                    continue;
                }

                // Define amino acid positions to show.

                let show_aa = build_show_aa(
                    ctl,
                    &rsi,
                    &vars_amino,
                    &shares_amino,
                    refdata,
                    dref,
                    &exacts,
                    exact_clonotypes,
                );

                // Define field types corresponding to the amino acid positions to show.

                let field_types = compute_field_types(ctl, &rsi, &show_aa);

                // Build varmat.

                let mut varmat = vec![vec![vec![b'-']; cols]; nexacts];
                for (col, (mat, (seqss, vars))) in mat
                    .iter()
                    .zip(rsi.seqss.iter().zip(vars.iter()))
                    .take(cols)
                    .enumerate()
                {
                    for (u, (m, seq)) in varmat.iter_mut().zip(mat.iter().zip(seqss.iter())) {
                        if m.is_some() {
                            let mut v = Vec::<u8>::new();
                            for &p in vars {
                                v.push(seq[p]);
                            }
                            u[col] = v;
                        }
                    }
                }

                // Find the fields associated to nd<k> if used.

                let mut lvarsc = lvars.clone();
                let mut nd_fields = Vec::<String>::new();
                for (i, x) in lvars.iter().enumerate() {
                    if x.starts_with("nd")
                        && x.after("nd").parse::<usize>().is_ok()
                        && x.after("nd").force_usize() >= 1
                    {
                        lvarsc.clear();
                        lvarsc.extend(lvars.iter().take(i).cloned());
                        let k = x.after("nd").force_usize();
                        let mut n = vec![0_usize; ctl.origin_info.n()];
                        for u in 0..nexacts {
                            let ex = &exact_clonotypes[exacts[u]];
                            for l in 0..ex.ncells() {
                                n[ex.clones[l][0].dataset_index] += 1;
                            }
                        }
                        let mut datasets = ctl.origin_info.dataset_id.clone();
                        // does not work for unknown reason, so "manually" replaced
                        // sort_sync2(&mut n, &mut datasets);
                        let permutation = permutation::sort(&n[..]);
                        n = permutation.apply_slice(&n[..]);
                        datasets = permutation.apply_slice(&datasets[..]);
                        n.reverse();
                        datasets.reverse();
                        for l in 0..n.len() {
                            if n[l] == 0 {
                                n.truncate(l);
                                datasets.truncate(l);
                                break;
                            }
                        }
                        for (l, ds) in datasets.iter().take(k).enumerate() {
                            if l >= n.len() {
                                break;
                            }
                            nd_fields.push(format!("n_{}", ds.as_str()));
                            lvarsc.push(format!("n_{}", ds.as_str()));
                        }
                        if n.len() > k {
                            nd_fields.push("n_other".to_string());
                            lvarsc.push("n_other".to_string());
                        }
                        lvarsc.extend(lvars.iter().skip(i + 1).cloned());
                        break;
                    }
                }
                let lvars = lvarsc.clone();
                let mut lvarsh = HashSet::<String>::new();
                for x in lvars.iter() {
                    lvarsh.insert(x.to_string());
                }

                // Now build table content.

                let mut sr = Vec::<(Vec<String>, Vec<Vec<String>>, Vec<Vec<u8>>, usize)>::new();
                let mut groups = HashMap::<usize, Vec<usize>>::new();
                for lvar in &lvars {
                    if let Some(Ok(d)) = lvar.strip_prefix('g').map(str::parse::<usize>) {
                        if groups.contains_key(&d) {
                            continue;
                        }
                        let mut e: EquivRel = EquivRel::new(nexacts as i32);
                        for (u1, &e1) in exacts.iter().take(nexacts).enumerate() {
                            let ex1 = &exact_clonotypes[e1];
                            for (u2, &e2) in exacts.iter().enumerate().take(nexacts).skip(u1 + 1) {
                                if e.class_id(u1 as i32) == e.class_id(u2 as i32) {
                                    continue;
                                }
                                let ex2 = &exact_clonotypes[e2];
                                let mut diffs = 0;
                                'comp: for (mm, vars) in mat.iter().zip(vars.iter()).take(cols) {
                                    if let (Some(m1), Some(m2)) = (mm[u1], mm[u2]) {
                                        let (s1, s2) =
                                            (&ex1.share[m1].seq_del, &ex2.share[m2].seq_del);
                                        for &p in vars {
                                            if s1[p] != s2[p] {
                                                diffs += 1;
                                                if diffs > d {
                                                    break 'comp;
                                                }
                                            }
                                        }
                                    }
                                }
                                if diffs <= d {
                                    e.join(u1 as i32, u2 as i32);
                                }
                            }
                        }
                        let mut c = Vec::<usize>::new();
                        let mut reps = Vec::<i32>::new();
                        e.orbit_reps(&mut reps);
                        for u in 0..nexacts {
                            c.push(bin_position(&reps, &e.class_id(u as i32)) as usize);
                        }
                        groups.insert(d, c);
                    }
                }

                // Set up to record stats that assign a value to each cell for a given variable.

                let mut stats = Vec::<(String, Vec<String>)>::new();

                // Compute some stats;

                let mut cred = Vec::<Vec<String>>::new();
                let mut pe = Vec::<Vec<String>>::new();
                let mut ppe = Vec::<Vec<String>>::new();
                let mut npe = Vec::<Vec<String>>::new();
                compute_some_stats(
                    ctl,
                    &lvars,
                    &exacts,
                    exact_clonotypes,
                    gex_info,
                    vdj_cells,
                    &n_vdj_gex,
                    &mut cred,
                    &mut pe,
                    &mut ppe,
                    &mut npe,
                );

                // Precompute for near and far.

                let mut fp = vec![Vec::<usize>::new(); varmat.len()]; // footprints
                for i in 0..varmat.len() {
                    for j in 0..varmat[i].len() {
                        if varmat[i][j] != vec![b'-'] {
                            fp[i].push(j);
                        }
                    }
                }

                // Form CDR3 consensus sequences.

                let mut cdr3_con = Vec::<Vec<u8>>::new();
                if ctl.gen_opt.color == "codon-diffs" {
                    cdr3_con = consensus_codon_cdr3(&rsi, &exacts, exact_clonotypes);
                }

                // Build rows.

                let mut cell_count = 0;
                for u in 0..nexacts {
                    let mut typex = vec![false; cols];
                    let mut row = Vec::<String>::new();
                    let mut gex_low = 0;
                    let mut cx = Vec::<Vec<String>>::new();
                    for col in 0..cols {
                        cx.push(vec![String::new(); rsi.cvars[col].len()]);
                    }
                    let clonotype_id = exacts[u];
                    let ex = &exact_clonotypes[clonotype_id];
                    let mut d_all = vec![Vec::<u32>::new(); ex.clones.len()];
                    let mut ind_all = vec![Vec::<u32>::new(); ex.clones.len()];
                    let mut these_stats = Vec::<(String, Vec<String>)>::new();
                    let resx = row_fill(
                        pass,
                        u,
                        ctl,
                        &exacts,
                        &mults,
                        exact_clonotypes,
                        gex_info,
                        refdata,
                        &varmat,
                        &fp,
                        &vars_amino,
                        &show_aa,
                        &ref_diff_pos,
                        &field_types,
                        &mut bads,
                        &mut gex_low,
                        &mut row,
                        &mut out_data,
                        &mut cx,
                        &mut d_all,
                        &mut ind_all,
                        &rsi,
                        dref,
                        &groups,
                        d_readers,
                        ind_readers,
                        h5_data,
                        &mut these_stats,
                        &stats_pass1,
                        vdj_cells,
                        &n_vdj_gex,
                        &lvars,
                        &lvarsh,
                        &nd_fields,
                        &peer_groups,
                        &extra_args,
                        &all_vars,
                        need_gex,
                        fate,
                        &cdr3_con,
                        allele_data,
                    );
                    stats.append(&mut these_stats.clone());
                    if pass == 1 {
                        stats_pass1.push(these_stats.clone());
                    }
                    these_stats.sort_by(|a, b| a.0.cmp(&b.0));
                    if let Err(e) = resx {
                        res.13 = e;
                        return;
                    }
                    let mut bli = ex
                        .clones
                        .iter()
                        .enumerate()
                        .map(|(l, clone)| (clone[0].barcode.clone(), clone[0].dataset_index, l))
                        .collect::<Vec<_>>();
                    // WHY ARE WE SORTING HERE?
                    bli.sort();
                    for col in 0..cols {
                        if mat[col][u].is_some() {
                            typex[col] = true;
                        }
                    }
                    for mut cxr in cx {
                        row.append(&mut cxr);
                    }
                    res.5 = gex_low;

                    // Compute per-cell entries.

                    if pass == 2 {
                        let mut subrows = Vec::<Vec<String>>::new();
                        compute_bu(
                            u,
                            cell_count,
                            &exacts,
                            &lvars,
                            ctl,
                            &bli,
                            ex,
                            exact_clonotypes,
                            &mut row,
                            &mut subrows,
                            &varmat,
                            have_gex,
                            gex_info,
                            &rsi,
                            &mut sr,
                            fate,
                            &nd_fields,
                            &alt_bcs,
                            &cred,
                            &pe,
                            &ppe,
                            &npe,
                            &d_all,
                            &ind_all,
                            mat,
                            &these_stats,
                            refdata,
                        );
                    }
                    cell_count += ex.clones.len();
                }
                let mut rord = Vec::<usize>::new(); // note that this is now superfluous
                for j in 0..sr.len() {
                    rord.push(j);
                }

                // Combine stats for the same variable.  This is needed because each exact
                // subclonotype contributes.  Note that we don't care about the order of the
                // values here (other than stability) because what we're going to do with them is
                // compute the mean or max.

                stats.sort_by(|a, b| a.0.cmp(&b.0));
                let stats_orig = stats.clone();
                let mut stats2 = Vec::<(String, Vec<String>)>::new();
                let mut i = 0;
                while i < stats.len() {
                    let mut j = i + 1;
                    while j < stats.len() {
                        if stats[j].0 != stats[i].0 {
                            break;
                        }
                        j += 1;
                    }
                    let all = stats[i..j]
                        .iter()
                        .flat_map(|s| s.1.iter().cloned())
                        .collect();
                    stats2.push((stats[i].0.clone(), all));
                    i = j;
                }
                stats = stats2;

                // Traverse the bounds and apply them.
                // Notes:
                // 1. This seems to run during both pass 1 and 2, and should only run
                //    during pass 1.
                // 2. The results of this can be counterintuitive, because the filtering is
                //    applied during pass 1, when there could be cells in the clonotype, that
                //    are removed by other filters.

                for bi in 0..ctl.clono_filt_opt.bounds.len() {
                    let x = &ctl.clono_filt_opt.bounds[bi];
                    let mut means = Vec::<f64>::new();
                    let mut mins = Vec::<f64>::new();
                    let mut maxs = Vec::<f64>::new();
                    // traverse the coefficients on the left hand side (each having a variable)
                    let mut fail = false;
                    for vi in x.var.iter().take(x.n()) {
                        let mut vals = Vec::<f64>::new(); // the stats for the variable
                        for stat in &stats {
                            if stat.0 == *vi {
                                for sk in &stat.1 {
                                    if let Ok(sk) = sk.parse::<f64>() {
                                        vals.push(sk);
                                    }
                                }
                                break;
                            }
                        }
                        let mut min = 1_000_000_000.0_f64;
                        let mut mean = 0.0;
                        let mut max = -1_000_000_000.0_f64;
                        let mut count = 0;
                        for val in vals {
                            if !val.is_nan() {
                                min = min.min(val);
                                mean += val;
                                max = max.max(val);
                                count += 1;
                            }
                        }
                        if count == 0 {
                            fail = true;
                        } else {
                            mins.push(min);
                            mean /= count as f64;
                            means.push(mean);
                            maxs.push(max);
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "mean" && (fail || !x.satisfied(&means))
                    {
                        if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                            in_center = false;
                        } else {
                            for b in bads.iter_mut().take(nexacts) {
                                *b = true;
                            }
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "min" && (fail || !x.satisfied(&mins)) {
                        if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                            in_center = false;
                        } else {
                            for b in bads.iter_mut().take(nexacts) {
                                *b = true;
                            }
                        }
                    }
                    if ctl.clono_filt_opt.bound_type[bi] == "max" && (fail || !x.satisfied(&maxs)) {
                        if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                            in_center = false;
                        } else {
                            for b in bads.iter_mut().take(nexacts) {
                                *b = true;
                            }
                        }
                    }
                }

                // Process COMPLETE.

                process_complete(ctl, nexacts, &mut bads, mat);

                // Done unless on second pass.

                if pass == 1 {
                    continue;
                }

                // See if we're in the test and control sets for gene scan.

                gene_scan_test(
                    ctl,
                    &stats,
                    &stats_orig,
                    nexacts,
                    n,
                    &mut res.9,
                    &mut res.10,
                );

                // Make the table.

                let mut logz = String::new();
                finish_table(
                    n,
                    ctl,
                    &exacts,
                    exact_clonotypes,
                    &rsi,
                    &vars,
                    &show_aa,
                    &field_types,
                    &lvars,
                    refdata,
                    dref,
                    &peer_groups,
                    &mut mlog,
                    &mut logz,
                    &stats,
                    &mut sr,
                    &extra_args,
                    pcols_sort,
                    &mut out_data,
                    &rord,
                    pass,
                    &cdr3_con,
                );

                // Save.

                res.1.push(logz);
                res.2.push((exacts.clone(), rsi.clone()));
                res.12.push(in_center);
                for u in 0..exacts.len() {
                    res.8 += exact_clonotypes[exacts[u]].ncells() as isize;
                }
            }
            if pass == 2 {
                res.7.push(out_data);
            }
        }
    });
    let exacts = exacts;
    for r in &results {
        if !r.13.is_empty() {
            return Err(r.13.clone());
        }
    }

    for ri in &results {
        for vj in &ri.11 {
            fate[vj.0].insert(vj.1.clone(), vj.2);
        }
    }

    // Sort results in descending order by number of cells.

    results.sort_by_key(|x| -x.8);

    // Write loupe output.

    let mut all_loupe_clonotypes = Vec::<Clonotype>::new();
    for r in results.iter_mut() {
        all_loupe_clonotypes.append(&mut r.6);
    }
    loupe_out(ctl, all_loupe_clonotypes, refdata, dref);

    // Write out the fate of each filtered cell. It would be better to collect
    // details of filtering in an enum as opposed to a string
    if !ctl.gen_opt.fate_file.is_empty() {
        let mut wtr = BufWriter::new(
            File::create(&ctl.gen_opt.fate_file).expect("Unable to open FATE_FILE for writing"),
        );
        serde_json::to_writer_pretty(&mut wtr, fate).map_err(|e| e.to_string())?;
    }

    // Set up to group and print clonotypes.

    for ri in results.iter_mut().take(orbits.len()) {
        for (v1, (v2, &v12)) in ri.1.iter().zip(ri.2.iter().zip(ri.12.iter())) {
            pics.push(v1.clone());
            exacts.push(v2.0.clone());
            rsi.push(v2.1.clone());
            in_center.push(v12);
        }
        out_datas.append(&mut ri.7);
    }

    // Gather some data for gene scan.

    if ctl.gen_opt.gene_scan_test.is_some() && !ctl.gen_opt.gene_scan_exact {
        for (i, r) in results.iter().take(orbits.len()).enumerate() {
            for (&v9, &v10) in r.9.iter().zip(r.10.iter()) {
                if v9 {
                    tests.push(i);
                }
                if v10 {
                    controls.push(i);
                }
            }
        }
    }
    if ctl.gen_opt.gene_scan_test.is_some() && ctl.gen_opt.gene_scan_exact {
        for (r, e) in results.iter().zip(exacts.iter()) {
            for (&ej, (&v9, &v10)) in e.iter().zip(r.9.iter().zip(r.10.iter())) {
                if v9 {
                    tests.push(ej);
                }
                if v10 {
                    controls.push(ej);
                }
            }
        }
    }
    Ok(())
}
