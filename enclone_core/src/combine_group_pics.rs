// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::defs::justification;
use ansi_escape::{emit_eight_bit_color_escape, emit_end_escape};
use io_utils::{fwrite, fwriteln};
use std::io::Write;
use string_utils::{stringme, strme};
use tables::print_tabular;

pub fn combine_group_pics(
    group_pics: &[String],
    last_widths: &[u32],
    parseable_stdouth: bool,
    noprint: bool,
    noprintx: bool,
    html: bool,
    ngroup: bool,
    pretty: bool,
) -> String {
    let mut glog = Vec::<u8>::new();
    let mut done = false;
    if noprint && parseable_stdouth && !group_pics.is_empty() {
        let mut rows = Vec::<Vec<String>>::new();
        for pic in group_pics {
            let r: Vec<_> = pic.split('\n').collect();
            for rj in r.iter().take(r.len() - 1) {
                let s = rj.split('\t').map(str::to_owned).collect();
                rows.push(s);
            }
        }
        let mut same = true;
        let n = rows[0].len();
        for row in rows.iter().skip(1) {
            if row.len() != n {
                same = false;
            }
        }
        if same {
            let justify = rows[0]
                .iter()
                .map(String::as_str)
                .map(justification)
                .collect();
            print_tabular(&mut glog, &rows, 2, Some(justify));
            done = true;
        }
    }

    if !done {
        // Get the newlines right is tricky, so they're marked.

        for i in 0..group_pics.len() {
            if !noprint {
                if !html && !ngroup && (!noprintx || i > 0) {
                    fwriteln!(glog, ""); // NEWLINE 1
                }

                // If we just printed a clonotype box, output a bar.

                if i > 0 && last_widths[i - 1] > 0 {
                    if ngroup || html {
                        fwriteln!(glog, ""); // NEWLINE 2
                    }
                    if pretty {
                        let mut log = Vec::<u8>::new();
                        emit_eight_bit_color_escape(&mut log, 44);
                        fwrite!(glog, "{}", strme(&log));
                    }
                    fwrite!(glog, "╺{}╸", "━".repeat((last_widths[i - 1] - 2) as usize));
                    if !ngroup {
                        fwriteln!(glog, ""); // NEWLINE 3
                    }
                    fwriteln!(glog, ""); // NEWLINE 4
                    if pretty {
                        let mut log = Vec::<u8>::new();
                        emit_end_escape(&mut log);
                        fwrite!(glog, "{}", strme(&log));
                    }
                }
            }
            glog.append(&mut group_pics[i].as_bytes().to_vec());
        }
    }
    stringme(&glog)
}
