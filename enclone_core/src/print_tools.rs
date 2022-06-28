// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use ansi_escape::{emit_end_escape, print_color};
use io_utils::fwrite;
use std::io::Write;
use string_utils::strme;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn emit_codon_color_escape(c: &[u8], log: &mut Vec<u8>) {
    let mut s = 0;
    if c == b"CTG" {
        s = 3;
    } else if c == b"AGG" {
        s = 1;
    } else if c == b"AGT" {
        s = 2;
    } else {
        for i in 0..3 {
            if c[i] == b'A' {
            } else if c[i] == b'C' {
                s += 1;
            } else if c[i] == b'G' {
                s += 2;
            } else if c[i] == b'T' {
                s += 3;
            } else {
                panic!("Illegal codon: \"{}\".", strme(c));
            }
        }
    }
    let s = s % 6;
    print_color(s, log);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn color_by_property(c: &[u8], log: &mut Vec<u8>) {
    for &ci in c {
        let mut color = 7;
        if ci == b'A' || ci == b'G' || ci == b'I' || ci == b'L' || ci == b'P' || ci == b'V' {
            color = 0;
        } else if ci == b'F' || ci == b'W' || ci == b'Y' {
            color = 1;
        } else if ci == b'D' || ci == b'E' {
            color = 2;
        } else if ci == b'R' || ci == b'H' || ci == b'K' {
            color = 3;
        } else if ci == b'S' || ci == b'T' {
            color = 4;
        } else if ci == b'C' || ci == b'M' {
            color = 5;
        } else if ci == b'N' || ci == b'Q' {
            color = 6;
        }
        if color < 7 {
            print_color(color, log);
        }
        fwrite!(log, "{}", ci as char);
        if color < 7 {
            emit_end_escape(log);
        }
    }
}
