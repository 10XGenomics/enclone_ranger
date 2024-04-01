// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::{EncloneControl, GexInfo};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn get_gex_matrix_entry(
    ctl: &EncloneControl,
    gex_info: &GexInfo,
    fid: usize,
    d_all: &[Vec<u32>],
    ind_all: &[Vec<u32>],
    li: usize,
    l: usize,
    y: &str,
) -> f64 {
    let mut raw_count = 0 as f64;

    for j in 0..d_all[l].len() {
        if ind_all[l][j] == fid as u32 {
            raw_count = d_all[l][j] as f64;
            break;
        }
    }

    let mult = if y.ends_with("_g") {
        gex_info.gex_mults[li]
    } else {
        gex_info.fb_mults[li]
    };
    if !ctl.gen_opt.full_counts {
        raw_count *= mult;
    }
    raw_count
}
