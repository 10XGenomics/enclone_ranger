// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file supplies the single function print_clonotypes.  It prints clonotypes, but also
// does some filtering to remove 'noise' clonotypes.
//
// Problem: stack traces from this file consistently do not go back to the main program.

use crate::define_mat::{define_mat, Od};
use crate::filter::survives_filter;
use crate::loupe::{loupe_out, make_loupe_clonotype};
use crate::print_utils3::define_column_info;
use crate::print_utils5::delete_weaks;
use enclone_core::defs::{CloneInfo, ColInfo};
use enclone_core::enclone_structs::{BarcodeFates, EncloneExacts, EncloneSetup, GexReaders};
use enclone_proto::types::Clonotype;
use itertools::Itertools;
use rayon::prelude::*;
use std::cmp::Reverse;
use std::fs::File;
use std::io::BufWriter;
use vector_utils::{erase_if, next_diff12_3};

/// Print clonotypes.  A key challenge here is to define the columns that represent shared
/// chains.  This is given below by the code that forms an equivalence relation on the CDR3_AAs.
///
/// This code carries out a second function, which is to filter out exact subclonotypes in orbits
/// that appear to be junk.  Exactly how these should be reflected in files is TBD.
///
/// Some inputs for this section:
/// refdata                = reference sequence data
/// ctl                    = control parameters
/// exact_clonotypes       = the vector of all exact subclonotypes
/// info                   = vector of clonotype info
/// eq                     = equivalence relation on info
pub fn print_clonotypes<T: Send>(
    setup: &EncloneSetup,
    enclone_exacts: &EncloneExacts,
    gex_readers: &[Option<GexReaders<'_>>],
    fate: &[BarcodeFates],
    mut proc: impl OrbitProcessor<T> + Send + Sync,
    // proc: impl OrbitProcessor<T>,
) -> Result<(), String> {
    let EncloneSetup {
        ctl,
        ann: _,
        gex_info,
        tall: _,
        refdata,
    } = setup;
    let EncloneExacts {
        to_bc,
        exact_clonotypes,
        raw_joins,
        info,
        orbits,
        vdj_cells: _,
        join_info: _,
        drefs: dref,
        sr,
        allele_data: _,
    } = enclone_exacts;

    // Traverse the orbits.

    // 0: index in reps
    // 1: vector of clonotype pictures
    // 2: vector of some clonotype info
    //    [parallel to 1]
    // next to last three entries = whitelist contam, denominator for that, low gex count
    // added out_datas (used to be next to last three, now one more)
    let result_iter = orbits.par_iter().map(|o| {
        let od: Vec<_> = o
            .iter()
            .map(|id| {
                let x: &CloneInfo = &info[*id as usize];
                (x.origin.clone(), x.clonotype_id, *id)
            })
            .sorted()
            .collect();

        // Reconstruct the participating clones.  This is needed because most exact subclonotypes
        // having more than two chains have been split up.
        //
        // Capture these data into parallel data structures, one per exact subclonotype:
        // exacts: the exact subclonotype ids
        // mults:  number of cells [redundant, might remove]

        let mut exacts = Vec::<usize>::new();
        let mut mults = Vec::<usize>::new();
        let mut j = 0;
        while j < od.len() {
            let k = next_diff12_3(&od, j);
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

        sort_exact_clonotypes(setup, enclone_exacts, &od, &mut exacts, &mut mults);

        // Define a matrix mat[col][ex] which is the column of the exact subclonotype
        // corresponding to the given column col of the clonotype, which may or may not be
        // defined.  Then define other information associated to each chain.  These are
        // reference sequence identifiers, CDR3 start positions, and the like.

        let mat = define_mat(
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
        let mut rsi = define_column_info(ctl, &exacts, exact_clonotypes, &mat, refdata);
        rsi.mat = mat;
        let mat = &rsi.mat;

        // Filter.

        // Let n be the total number of cells in this pass.
        let n: usize = mults.iter().sum();

        if n >= ctl.clono_filt_opt.ncells_low
            || ctl.clono_group_opt.asymmetric_center == "from_filters"
        {
            // Mark some weak exact subclonotypes for deletion.
            delete_weaks(ctl, &exacts, exact_clonotypes, mat, refdata, &mut bads);

            proc.filter(
                setup,
                enclone_exacts,
                gex_readers,
                fate,
                &exacts,
                &mults,
                n,
                &rsi,
                &mut bads,
                &mut stats_pass1,
                true,
            )?;
        }

        // Delete weak exact subclonotypes.

        if !ctl.clono_filt_opt.protect_bads {
            erase_if(&mut mults, &bads);
            erase_if(&mut exacts, &bads);
        }

        sort_exact_clonotypes(setup, enclone_exacts, &od, &mut exacts, &mut mults);

        // Define a matrix mat[col][ex] which is the column of the exact subclonotype
        // corresponding to the given column col of the clonotype, which may or may not be
        // defined.  Then define other information associated to each chain.  These are
        // reference sequence identifiers, CDR3 start positions, and the like.

        let mat = define_mat(
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
        let mut rsi = define_column_info(ctl, &exacts, exact_clonotypes, &mat, refdata);
        rsi.mat = mat;

        // Filter.

        let mut in_center = true;
        if !survives_filter(
            &exacts,
            &rsi,
            ctl,
            exact_clonotypes,
            refdata,
            gex_info,
            dref,
        ) {
            if ctl.clono_group_opt.asymmetric_center == "from_filters" {
                in_center = false;
            } else {
                return Ok((0, None, None));
            }
        }

        // Generate Loupe data.

        let loupe_clonotype = (!ctl.gen_opt.binary.is_empty() || !ctl.gen_opt.proto.is_empty())
            .then(|| make_loupe_clonotype(exact_clonotypes, &exacts, &rsi, refdata, dref, ctl));

        // Let n be the total number of cells in this pass.

        let n: usize = mults.iter().sum();

        if !(n >= ctl.clono_filt_opt.ncells_low
            || ctl.clono_group_opt.asymmetric_center == "from_filters")
        {
            return Ok((0, loupe_clonotype, None));
        }

        let num_cells: usize = exacts
            .iter()
            .map(|exact| exact_clonotypes[*exact].ncells())
            .sum();

        let res = proc.finalize(
            setup,
            enclone_exacts,
            gex_readers,
            fate,
            &exacts,
            &mults,
            n,
            &rsi,
            &mut bads,
            &mut stats_pass1,
            in_center,
        )?;

        Ok((num_cells, loupe_clonotype, res))
    });
    let mut results: Vec<_> =
        result_iter.collect::<Result<Vec<(usize, Option<Clonotype>, Option<T>)>, String>>()?;

    // Sort results in descending order by number of cells.

    results.sort_by_key(|(num_cells, _, _)| Reverse(*num_cells));

    // Write out the fate of each filtered barcode.
    if !ctl.gen_opt.fate_file.is_empty() {
        let mut wtr = BufWriter::new(
            File::create(&ctl.gen_opt.fate_file).expect("Unable to open FATE_FILE for writing"),
        );
        serde_json::to_writer_pretty(&mut wtr, &fate).map_err(|e| e.to_string())?;
    }

    let mut all_loupe_clonotypes = Vec::<Clonotype>::new();

    for (_, loupe_clonotype, ri) in results {
        all_loupe_clonotypes.extend(loupe_clonotype);
        proc.collect(ri);
    }

    // Write loupe output.
    loupe_out(ctl, all_loupe_clonotypes, refdata, dref);

    Ok(())
}

/// Sort exact subclonotypes.
fn sort_exact_clonotypes(
    setup: &EncloneSetup,
    enclone_exacts: &EncloneExacts,
    od: &[Od],
    exacts: &mut Vec<usize>,
    mults: &mut Vec<usize>,
) {
    let EncloneSetup {
        ctl,
        ann: _,
        gex_info: _,
        tall: _,
        refdata,
    } = setup;
    let EncloneExacts {
        to_bc,
        exact_clonotypes,
        raw_joins,
        info,
        orbits: _,
        vdj_cells: _,
        join_info: _,
        drefs: dref,
        sr,
        allele_data: _,
    } = enclone_exacts;
    let mat = define_mat(
        to_bc,
        sr,
        ctl,
        exact_clonotypes,
        exacts,
        od,
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
    *exacts = permutation.apply_slice(&exacts[..]);
    *mults = permutation.apply_slice(&mults[..]);
    exacts.reverse();
    mults.reverse();
}

/// Inject a behavior to provide additional filtering and post-processing of each orbit.
pub trait OrbitProcessor<T> {
    #[allow(unused)]
    fn filter(
        &self,
        setup: &EncloneSetup,
        enclone_exacts: &EncloneExacts,
        gex_readers: &[Option<GexReaders<'_>>],
        fate: &[BarcodeFates],
        exacts: &[usize],
        mults: &[usize],
        n: usize,
        rsi: &ColInfo,
        bads: &mut [bool],
        stats_pass1: &mut Vec<Vec<(String, Vec<String>)>>,
        in_center: bool,
    ) -> Result<(), String> {
        Ok(())
    }

    #[allow(unused)]
    fn finalize(
        &self,
        setup: &EncloneSetup,
        enclone_exacts: &EncloneExacts,
        gex_readers: &[Option<GexReaders<'_>>],
        fate: &[BarcodeFates],
        exacts: &[usize],
        mults: &[usize],
        n: usize,
        rsi: &ColInfo,
        bads: &mut [bool],
        stats_pass1: &mut Vec<Vec<(String, Vec<String>)>>,
        in_center: bool,
    ) -> Result<Option<T>, String> {
        Ok(None)
    }

    /// Collect a single result from processing an orbit; the results will be
    /// provided in a sorted order.
    #[allow(unused)]
    fn collect(&mut self, result: Option<T>) {}
}
