// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Extract the entries in a given all_contig_annotations.json file that corrrespond to barcodes
// in a given sorted vector.

use io_utils::open_userfile_for_read;
use std::fmt::Write;
use std::io::BufRead;
use string_utils::TextUtils;
use vector_utils::bin_member;

pub fn subset_all_contig_annotations_json(filename: &str, barcodes: &[String]) -> String {
    let mut x = "[\n".to_string();
    let mut lines = Vec::<String>::new();
    let f = open_userfile_for_read(filename);
    let mut keep = false;
    for line in f.lines() {
        let s = line.unwrap();
        if s.starts_with('[') {
            continue;
        }
        if s == "]" {
            if keep {
                for line in &lines {
                    writeln!(x, "{line}").unwrap();
                }
            }
            break;
        }
        lines.push(s);
        let s = lines.last().unwrap().as_str();
        if s.starts_with("        \"barcode\": \"") {
            let t = s.between("        \"barcode\": \"", "\"");
            keep = bin_member(barcodes, &t.to_string());
        } else if s.starts_with("    }") {
            if keep {
                for line in &lines {
                    writeln!(&mut x, "{line}").unwrap();
                }
            }
            lines.clear();
        }
    }
    x += "]\n";
    x
}
