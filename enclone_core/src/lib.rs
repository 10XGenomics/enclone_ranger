// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod align_to_vdj_ref;
pub mod allowed_vars;
pub mod barcode_fate;
pub mod cell_color;
pub mod combine_group_pics;
pub mod defs;
pub mod enclone_structs;
pub mod h5_lazy;
pub mod hcomp;
pub mod join_one;
pub mod linear_condition;
pub mod logging;
pub mod median;
pub mod opt_d;
pub mod print_tools;
pub mod slurp;
pub mod stringulate;
pub mod test_def;
pub mod var_reg;

use io_utils::path_exists;
use std::cmp::max;
use std::fmt::Write;
use std::fs::{remove_file, File};
use std::io::BufRead;
use string_utils::TextUtils;

#[cfg(not(target_os = "windows"))]
use string_utils::stringme;

#[cfg(not(target_os = "windows"))]
use tilde_expand::tilde_expand;

// tilde_expand_me: not that this is NOT implementd for Windows

pub fn tilde_expand_me(s: &mut String) {
    #[cfg(not(target_os = "windows"))]
    {
        *s = stringme(&tilde_expand(s.as_bytes()));
    }
}

pub fn hcat(col1: &[String], col2: &[String], sep: usize) -> Vec<String> {
    let mut cat = Vec::<String>::new();
    let height = max(col1.len(), col2.len());
    let mut width1 = 0;
    for x in col1 {
        width1 = max(width1, x.chars().count() + sep);
    }
    for i in 0..height {
        let mut s = if i < col1.len() {
            col1[i].clone()
        } else {
            String::new()
        };
        while s.chars().count() < width1 {
            s += " ";
        }
        if i < col2.len() {
            s += &col2[i];
        }
        cat.push(s);
    }
    cat
}

pub fn expand_integer_ranges(x: &str) -> String {
    let mut tokens = Vec::<String>::new();
    let mut token = String::new();
    for c in x.chars() {
        if c == ',' || c == ':' || c == ';' {
            if !token.is_empty() {
                tokens.push(token.clone());
                token.clear();
            }
            tokens.push(c.to_string());
        } else {
            token.push(c);
        }
    }
    if !token.is_empty() {
        tokens.push(token);
    }
    let mut tokens2 = String::new();
    for token in tokens {
        if let Some((n1s, n2s)) = token.split_once('-') {
            if let Ok(n1) = n1s.parse::<usize>() {
                if let Ok(n2) = n2s.parse::<usize>() {
                    if n1 <= n2 {
                        for n in n1..=n2 {
                            if n > n1 {
                                tokens2.push(',');
                            }
                            write!(tokens2, "{n}").unwrap();
                        }
                        continue;
                    }
                }
            }
        }
        tokens2 += token.as_str();
    }
    tokens2
}

// Parse a line, breaking at blanks, but not if they're in quotes.  And strip the quotes.
// Ridiculously similar to parse_csv, probably should refactor.
pub fn parse_bsv(x: &str) -> Vec<&str> {
    let mut args = Vec::<&str>::new();
    let w = x.as_bytes();
    let (mut quotes, mut i) = (0, 0);
    while i < w.len() {
        let mut j = i;
        while j < w.len() {
            if quotes % 2 == 0 && w[j] == b' ' {
                break;
            }
            if w[j] == b'"' {
                quotes += 1;
            }
            j += 1;
        }
        // These will always be a char boundaries because it's either the end of
        // the string or it's a space character, which is a single byte in
        // UTF-8.  That remains true if i/j are '"' bytes.
        let (mut start, mut stop) = (i, j);
        if stop - start >= 2 && w[start] == b'"' && w[stop - 1] == b'"' {
            start += 1;
            stop -= 1;
        }
        args.push(&x[start..stop]);
        i = j + 1;
    }
    args
}

/// Test to see if a line can be read from the given file f.  If not, return an error message
/// the references arg, which is supposed to be the name of a command line argument from which
/// f originated.
pub fn require_readable_file(f: &str, arg: &str) -> Result<(), String> {
    let x = std::fs::File::open(f);
    if x.is_err() {
        return Err(format!(
            "\nThe file {} could not be opened because {}.\nThis came from \
            the command line argument {}.\n",
            f,
            x.err().unwrap(),
            arg,
        ));
    }
    let y = std::io::BufReader::new(x.unwrap());
    if let Some(line) = y.lines().next() {
        if line.is_err() {
            let mut err = line.err().unwrap().to_string();
            if err.starts_with("Is a directory") {
                err = "it is a directory".to_string();
            }
            return Err(format!(
                "\nThe file {f} could not be read because {err}.\nThis came from \
                the command line argument {arg}.\n",
            ));
        }
    }
    Ok(())
}

/// Test a file for writeability by writing and then deleting it.
pub fn test_writeable(val: &str, evil_eye: bool) -> Result<(), String> {
    if evil_eye {
        println!("creating file {val} to test writability");
    }
    let f = File::create(val);
    if f.is_err() {
        let mut msgx =
            format!("\nYou've specified an output file\n{val}\nthat cannot be written.\n");
        if val.contains('/') {
            let dir = val.rev_before("/");
            let msg = if path_exists(dir) {
                "exists"
            } else {
                "does not exist"
            };
            writeln!(msgx, "Note that the path {dir} {msg}.").unwrap();
        }
        return Err(msgx);
    }
    if evil_eye {
        println!("removing file {val}");
    }
    remove_file(val).unwrap_or_else(|_| panic!("could not remove file {val}"));
    if evil_eye {
        println!("removal of file {val} complete");
    }
    Ok(())
}
