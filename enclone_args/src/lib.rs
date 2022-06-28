// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.
#![allow(clippy::needless_range_loop)]

use io_utils::path_exists;

pub mod load_gex;
pub mod load_gex_core;
pub mod load_gex_util;
pub mod proc_args;
pub mod proc_args2;
pub mod proc_args3;
pub mod proc_args_check;
pub mod proc_args_post;
pub mod process_special_arg1;
pub mod process_special_arg2;
pub mod read_json;

// parse_csv_pure: same as parse_csv, but don't strip out quotes

pub fn parse_csv_pure(x: &str) -> Vec<&str> {
    let w = x.char_indices().collect::<Vec<_>>();
    let mut y = Vec::new();
    let (mut quotes, mut i) = (0, 0);
    while i < w.len() {
        let mut j = i;
        while j < w.len() {
            if quotes % 2 == 0 && w[j].1 == ',' {
                break;
            }
            if w[j].1 == '"' {
                quotes += 1;
            }
            j += 1;
        }
        let (start, stop) = (i, j);
        y.push(&x[w[start].0..w[stop].0]);
        i = j + 1;
    }
    if !w.is_empty() && w.last().unwrap().1 == ',' {
        y.push("");
    }
    y
}

pub fn fnx(outs: &str, name: &str) -> String {
    let mut file = format!("{}/../{}", outs, name);
    if !path_exists(&file) {
        file = format!("{}/{}", outs, name);
    }
    file
}
