// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file provides the tail end code for join.rs, plus a small function used there.

use enclone_core::defs::{CloneInfo, EncloneControl};
use equiv::EquivRel;
use stats_utils::percent_ratio;
use std::time::Instant;
use vector_utils::next_diff1_2;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn finish_join(
    ctl: &EncloneControl,
    info: &[CloneInfo],
    results: &[(
        usize,
        usize,
        usize,
        usize,
        Vec<(usize, usize, bool, Vec<u8>)>,
        Vec<(usize, usize)>,
    )],
    join_info: &mut Vec<(usize, usize, bool, Vec<u8>)>,
) -> EquivRel {
    // Tally results.

    let (mut joins, mut errors) = (0, 0);
    let timer3 = Instant::now();
    for r in results {
        joins += r.2;
        errors += r.3;
        for i in &r.4 {
            let u1 = i.0;
            let u2 = i.1;
            let err = i.2;
            let log = i.3.clone();
            join_info.push((u1, u2, err, log));
        }
    }
    if !ctl.silent {
        println!("{joins} joins");
        if ctl.origin_info.donors > 1 {
            println!("{errors} errors");
        }
    }

    // Make equivalence relation.

    let mut eq: EquivRel = EquivRel::new(info.len() as i32);
    for r in results {
        for j in &r.5 {
            eq.join(j.0 as i32, j.1 as i32);
        }
    }

    // Join orbits that cross subclones of a clone.  This arose because we split up multi-chain
    // clonotypes into two-chain clonotypes.

    let mut ox = Vec::<(usize, i32)>::new();
    for (i, x) in info.iter().enumerate() {
        ox.push((x.clonotype_id, eq.class_id(i as i32)));
    }
    ox.sort_unstable();
    let mut i = 0;
    while i < ox.len() {
        let j = next_diff1_2(&ox, i as i32) as usize;
        for k in i..j - 1 {
            eq.join(ox[k].1, ox[k + 1].1);
        }
        i = j;
    }

    // Tally whitelist contamination.
    // WARNING: THIS ONLY WORKS IF YOU RUN WITH CLONES=1 AND NO OTHER FILTERS.

    let mut white = ctl.clono_filt_opt_def.whitef;
    for j in 0..ctl.clono_print_opt.cvars.len() {
        if ctl.clono_print_opt.cvars[j] == "white" {
            white = true;
        }
    }
    if white {
        let mut bads = 0;
        let mut denom = 0;
        for r in results {
            bads += r.2;
            denom += r.3;
        }
        let bad_rate = percent_ratio(bads, denom);
        println!("whitelist contamination rate = {bad_rate:.2}%");
    }
    ctl.perf_stats(&timer3, "in tail of join");
    eq
}
