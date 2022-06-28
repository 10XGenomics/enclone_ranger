// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

pub mod align_to_vdj_ref;
pub mod allowed_vars;
pub mod cell_color;
pub mod combine_group_pics;
pub mod defs;
pub mod enclone_structs;
pub mod hcomp;
pub mod join_one;
pub mod linear_condition;
pub mod logging;
pub mod main_testlist;
pub mod mammalian_fixed_len;
pub mod median;
pub mod opt_d;
pub mod packing;
pub mod prepare_for_apocalypse;
pub mod print_tools;
pub mod set_speakers;
pub mod slurp;
pub mod stringulate;
pub mod test_def;
pub mod var_reg;

use lazy_static::lazy_static;
use std::cmp::max;
use std::env;
use std::fmt::Write;
use std::io::BufRead;
use std::sync::Mutex;
use std::time::Duration;

#[cfg(not(target_os = "windows"))]
use string_utils::stringme;

#[cfg(not(target_os = "windows"))]
use tilde_expand::tilde_expand;

// tilde_expand_me: not that this is NOT implementd for Windows

pub fn tilde_expand_me(_s: &mut String) {
    #[cfg(not(target_os = "windows"))]
    {
        *_s = stringme(&tilde_expand(&*_s.as_bytes()));
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
                            write!(tokens2, "{}", n).unwrap();
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

lazy_static! {
    pub static ref BUG_REPORT_ADDRESS: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
    pub static ref REMOTE_HOST: Mutex<Vec<String>> = Mutex::new(Vec::<String>::new());
}

const VERSION_STRING: &str = env!("VERSION_STRING");

// WARNING: the version string will not be up to date unless enclone_core/build.rs is touched
// before compiling.  Use current_version_string() to always get the current version.

pub fn version_string() -> String {
    VERSION_STRING.to_string()
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

pub fn fetch_url(url: &str) -> Result<String, String> {
    const TIMEOUT: u64 = 120; // timeout in seconds
    let req = attohttpc::get(url).read_timeout(Duration::new(TIMEOUT, 0));
    let response = req.send();
    if response.is_err() {
        return Err(format!(
            "\nFailed to access URL {},\ntimeout after two minutes.  There are a few ways that \
            you might have arrived at this state:\n• The server for that URL is down.\n\
            • The server for that URL is overloaded and responding very slowly.\n\
            • Same thing as last, and your process is slamming the server.  Please inspect \
            your command!\n\
            • There is a bug in this program.  This is relatively unlikely but possible.\n",
            url
        ));
    }
    let response = response.unwrap();
    if !response.is_success() {
        let msg = response.text().unwrap();
        if msg.contains("Not found") {
            return Err(format!(
                "\nAttempt to access the URL\n{}\nfailed with \"Not found\".  Could there \
                be something wrong with the id?\n",
                url
            ));
        }
        return Err(format!("Failed to access URL {}: {}.", url, msg));
    }
    Ok(response.text().unwrap())
}

// Test to see if a line can be read from the given file f.  If not, return an error message
// the references arg, which is supposed to be the name of a command line argument from which
// f originated.

pub fn require_readable_file(f: &str, arg: &str) -> Result<(), String> {
    let x = std::fs::File::open(&f);
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
                "\nThe file {} could not be read because {}.\nThis came from \
                the command line argument {}.\n",
                f, err, arg,
            ));
        }
    }
    Ok(())
}
