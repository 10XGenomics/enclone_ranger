// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::GeneScanOpts;

pub fn gene_scan_test(
    opts: &GeneScanOpts,
    exact: bool,
    stats: &[(String, Vec<String>)],
    stats_orig: &[(String, Vec<String>)],
    nexacts: usize,
    n: usize,
    in_test: &mut Vec<bool>,
    in_control: &mut Vec<bool>,
) {
    // See if we're in the test and control sets for gene scan (non-exact case).

    if !exact {
        let x = &opts.test;
        let means = x
            .var
            .iter()
            .take(x.n())
            .map(|xn| {
                stats
                    .iter()
                    .find_map(|stat| {
                        if stat.0 == *xn {
                            Some(
                                stat.1
                                    .iter()
                                    .filter_map(|k| k.parse::<f64>().ok())
                                    .sum::<f64>(),
                            )
                        } else {
                            None
                        }
                    })
                    .unwrap_or_default()
                    / n as f64
            })
            .collect::<Vec<_>>();

        in_test.push(x.satisfied(&means));
        let x = &opts.control;
        let means = x
            .var
            .iter()
            .take(x.n())
            .map(|xn| {
                stats
                    .iter()
                    .find_map(|stat| {
                        if stat.0 == *xn {
                            Some(
                                stat.1
                                    .iter()
                                    .filter_map(|k| k.parse::<f64>().ok())
                                    .sum::<f64>(),
                            )
                        } else {
                            None
                        }
                    })
                    .unwrap_or_default()
                    / n as f64
            })
            .collect::<Vec<_>>();
        in_control.push(x.satisfied(&means));
    } else {
        // See if we're in the test and control sets for gene scan (exact case).
        let x = &opts.test;
        for k in 0..nexacts {
            let mut means = Vec::<f64>::new();
            for xn in x.var.iter().take(x.n()) {
                let mut vals = Vec::<f64>::new();
                let mut count = 0;
                for stat in stats_orig {
                    if stat.0 == *xn {
                        if count == k {
                            for k in &stat.1 {
                                if let Ok(v) = k.parse::<f64>() {
                                    vals.push(v);
                                }
                            }
                            break;
                        }
                        count += 1;
                    }
                }
                let n = vals.len() as f64;
                means.push(vals.into_iter().sum::<f64>() / n);
            }
            in_test.push(x.satisfied(&means));
            let x = &opts.control;
            let mut means = Vec::<f64>::new();
            for xn in x.var.iter().take(x.n()) {
                let mut vals = Vec::<f64>::new();
                let mut count = 0;
                for stat in stats_orig {
                    if stat.0 == *xn {
                        if count == k {
                            for k in &stat.1 {
                                if let Ok(v) = k.parse::<f64>() {
                                    vals.push(v);
                                }
                            }
                            break;
                        }
                        count += 1;
                    }
                }
                means.push(vals.into_iter().sum::<f64>() / n as f64);
            }
            in_control.push(x.satisfied(&means));
        }
    }
}
