// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

// Functions print_tabular and print_tabular_vbox for making pretty tables.  And related utilities.

use io_utils::eprintme;
use itertools::Itertools;
use std::cmp::{max, min};
use string_utils::strme;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Package characters with ANSI escape codes that come before them.

pub fn package_characters_with_escapes(c: &[u8]) -> Vec<Vec<u8>> {
    let mut x = Vec::<Vec<u8>>::new();
    let mut escaped = false;
    let mut package = Vec::<u8>::new();
    for b in c {
        if escaped && *b != b'm' {
            package.push(*b);
        } else if *b == b'' {
            escaped = true;
            package.push(*b);
        } else if escaped && *b == b'm' {
            escaped = false;
            package.push(*b);
        } else {
            package.push(*b);
            x.push(package.clone());
            package.clear();
        }
    }
    x
}

pub fn package_characters_with_escapes_char(c: &[char]) -> Vec<Vec<char>> {
    let mut x = Vec::<Vec<char>>::new();
    let mut escaped = false;
    let mut package = Vec::<char>::new();
    for b in c {
        if escaped && *b != 'm' {
            package.push(*b);
        } else if *b == '' {
            escaped = true;
            package.push(*b);
        } else if escaped && *b == 'm' {
            escaped = false;
            package.push(*b);
        } else {
            package.push(*b);
            x.push(package.clone());
            package.clear();
        }
    }
    x
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Print out a matrix, with left-justified entries, and given separation between
// columns.  (Justification may be changed by supplying an optional argument
// consisting of a string of l's and r's.)

pub fn print_tabular(
    log: &mut Vec<u8>,
    rows: &[Vec<String>],
    sep: usize,
    justify: Option<Vec<u8>>,
) {
    let just = match justify {
        Some(x) => x,
        None => Vec::<u8>::new(),
    };
    let nrows = rows.len();
    let mut ncols = 0;
    for row in rows.iter().take(nrows) {
        ncols = max(ncols, row.len());
    }
    let mut maxcol = vec![0; ncols];
    for row in rows {
        for (j, item) in row.iter().enumerate() {
            maxcol[j] = max(maxcol[j], item.chars().count());
        }
    }
    for row in rows {
        for (j, x) in row.iter().enumerate() {
            if j < just.len() && just[j] == b'r' {
                log.append(&mut vec![b' '; maxcol[j] - x.chars().count()]);
                log.append(&mut x.as_bytes().to_vec());
                if j < row.len() - 1 {
                    log.append(&mut vec![b' '; sep]);
                }
            } else {
                log.append(&mut x.as_bytes().to_vec());
                if j < row.len() - 1 {
                    log.append(&mut vec![b' '; maxcol[j] - x.chars().count() + sep]);
                }
            }
        }
        log.push(b'\n');
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Compute the visible length of a string, counting unicode characters as width one and
// ignoring some ASCII escape sequences.

pub fn visible_width(s: &str) -> usize {
    let mut n = 0;
    let mut escaped = false;
    for c in s.chars() {
        if escaped && c != 'm' {
        } else if c == '' {
            escaped = true;
        } else if escaped && c == 'm' {
            escaped = false;
        } else {
            n += 1;
        }
    }
    n
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Print out a matrix, with given separation between columns.  Rows of the matrix
// may contain arbitrary UTF-8 and some escape sequences.  Put the entire thing in a box, with
// extra vertical bars.  The argument justify consists of symbols l and r, denoting
// left and right justification for given columns, respectively, and the symbol | to
// denote a vertical bar.
//
// There is no separation printed on the far left or far right.
//
// By a "matrix entry", we mean one of the Strings in "rows".
//
// Entries that begin with a backslash are reserved for future features.
// Symbols other than l or r or | in "justify" are reserved for future features.
//
// An entry may be followed on the right by one more entries whose contents are
// exactly "\ext".  In that case the entries are treated as multi-column.  Padding
// is inserted as needed on the "right of the multicolumn".
//
// An entry may be "\hline", which gets you a horizontal line.  The normal use case is to
// use one or more of these in succession horizontally to connect two vertical lines.  Cannot
// be combined with \ext.
//
// bold_box: use bold box characters
//
// Really only guaranteed to work for the tested cases.

pub fn print_tabular_vbox(
    log: &mut String,
    rows: &[Vec<String>],
    sep: usize,
    justify: &[u8],
    debug_print: bool,
    bold_box: bool,
) {
    // Define box characters.

    let dash = if !bold_box { '─' } else { '━' };
    let verty = if !bold_box { '│' } else { '┃' };
    let topleft = if !bold_box { '┌' } else { '┏' };
    let topright = if !bold_box { '┐' } else { '┓' };
    let botleft = if !bold_box { '└' } else { '┗' };
    let botright = if !bold_box { '┘' } else { '┛' };
    let tee = if !bold_box { '┬' } else { '┳' };
    let uptee = if !bold_box { '┴' } else { '┻' };
    let cross = if !bold_box { '┼' } else { '╋' };
    let lefty = if !bold_box { '├' } else { '┣' };
    let righty = if !bold_box { '┤' } else { '┫' };

    // Proceed.

    let mut rrr = rows.to_owned();
    let nrows = rrr.len();
    let mut ncols = 0;
    for item in rrr.iter().take(nrows) {
        ncols = max(ncols, item.len());
    }
    let mut vert = vec![false; ncols];
    let mut just = Vec::<u8>::new();
    let mut count = 0_isize;
    for item in justify {
        if *item == b'|' {
            assert!(count > 0);
            if count >= ncols as isize {
                eprintln!("\nposition of | in justify string is illegal");
                eprintme!(count, ncols);
            }
            assert!(count < ncols as isize);
            vert[(count - 1) as usize] = true;
        } else {
            just.push(*item);
            count += 1;
        }
    }
    if just.len() != ncols {
        eprintln!(
            "\nError.  Your table has {} columns but the number of \
             l or r symbols in justify is {}.\nThese numbers should be equal.",
            ncols,
            just.len()
        );
        eprintln!("justify = {}", strme(justify));
        for (i, row) in rows.iter().enumerate() {
            eprintln!("row {} = {} = {}", i + 1, row.len(), row.iter().format(","));
        }
        assert_eq!(just.len(), ncols);
    }
    let mut maxcol = vec![0; ncols];
    let mut ext = vec![0; ncols];
    for row in &rrr {
        for (j, item) in row.iter().enumerate() {
            if j < row.len() - 1 && row[j + 1] == "\\ext" {
                continue;
            }
            if item == "\\ext" || item == "\\hline" {
                continue;
            }
            maxcol[j] = max(maxcol[j], visible_width(item));
        }
    }
    if debug_print {
        println!("maxcol = {}", maxcol.iter().format(","));
    }

    // Add space according to ext entries.

    for i in 0..rrr.len() {
        for j in 0..rrr[i].len() {
            if j < rrr[i].len() - 1 && rrr[i][j + 1] == *"\\ext" && rrr[i][j] != *"\\ext" {
                let mut k = j + 1;
                while k < rrr[i].len() {
                    if rrr[i][k] != *"\\ext" {
                        break;
                    }
                    k += 1;
                }
                let need = visible_width(&rrr[i][j]);
                let mut have = 0;
                for l in j..k {
                    have += maxcol[l];
                    if l < k - 1 {
                        have += sep;
                        if vert[l] {
                            have += sep + 1;
                        }
                    }
                }
                if debug_print {
                    println!("row {i} column {j}, have = {have}, need = {need}");
                }
                if have > need {
                    if debug_print {
                        println!(
                            "adding {} spaces to right of row {} col {}",
                            have - need,
                            i,
                            j
                        );
                    }
                    for _ in need..have {
                        rrr[i][j].push(' ');
                    }
                } else if need > have {
                    maxcol[k - 1] += need - have;
                    if debug_print {
                        println!("increasing maxcol[{}] to {}", k - 1, maxcol[k - 1]);
                    }
                    ext[k - 1] += need - have;
                }
                let mut m = 0;
                for (u, row) in rrr.iter().enumerate() {
                    if j >= row.len() {
                        eprintln!("\nProblem with line {u}, not enough fields.\n");
                    }
                    if row[j] != *"\\ext" {
                        m = max(m, visible_width(&rrr[u][j]));
                    }
                }
                if m > visible_width(&rrr[i][j]) {
                    for _ in visible_width(&rrr[i][j])..m {
                        rrr[i][j].push(' ');
                    }
                }
            }
        }
    }

    // Create top boundary of table.

    log.push(topleft);
    for i in 0..ncols {
        let mut n = maxcol[i];
        if i < ncols - 1 {
            n += sep;
        }
        for _ in 0..n {
            log.push(dash);
        }
        if vert[i] {
            log.push(tee);
            for _ in 0..sep {
                log.push(dash);
            }
        }
    }
    log.push(topright);
    log.push('\n');

    // Go through the rows.

    for (i, row) in rrr.iter().take(nrows).enumerate() {
        if debug_print {
            println!("now row {} = {}", i, row.iter().format(","));
            println!("0 - pushing │ onto row {i}");
        }
        log.push(verty);
        for j in 0..min(ncols, row.len()) {
            // Pad entries according to justification.

            let mut x = String::new();
            if j >= row.len() {
                for _ in 0..maxcol[j] {
                    x.push(' ');
                }
            } else if row[j] == *"\\hline" {
                for _ in 0..maxcol[j] {
                    x.push(dash);
                }
            } else {
                let r = row[j].clone();
                let rlen = visible_width(&r);
                let mut xlen = 0;
                if r != *"\\ext" {
                    if just[j] == b'r' {
                        for _ in rlen..(maxcol[j] - ext[j]) {
                            x.push(' ');
                            xlen += 1;
                        }
                    }
                    if j < row.len() {
                        x += &r;
                        xlen += visible_width(&r);
                    }
                    if just[j] == b'r' {
                        for _ in (maxcol[j] - ext[j])..maxcol[j] {
                            x.push(' ');
                            xlen += 1;
                        }
                    }
                    if just[j] == b'l' {
                        for _ in xlen..maxcol[j] {
                            x.push(' ');
                        }
                    }
                }
            }
            for c in x.chars() {
                log.push(c);
            }

            // Add separations and separators.

            let mut add_sep = true;
            if j + 1 < row.len() && row[j + 1] == *"\\ext" {
                add_sep = false;
            }
            let mut jp = j;
            while jp + 1 < row.len() {
                if row[jp + 1] != *"\\ext" {
                    break;
                }
                jp += 1;
            }
            if add_sep && jp < ncols - 1 {
                if row[j] == *"\\hline" {
                    for _ in 0..sep {
                        log.push(dash);
                    }
                } else {
                    for _ in 0..sep {
                        log.push(' ');
                    }
                }
            }
            if vert[j] && row[j + 1] != "\\ext" {
                if debug_print {
                    println!("1 - pushing {verty} onto row {i}, j = {j}");
                }
                log.push(verty);
                if row[j + 1] == *"\\hline" {
                    for _ in 0..sep {
                        log.push(dash);
                    }
                } else {
                    for _ in 0..sep {
                        log.push(' ');
                    }
                }
            }
        }
        if debug_print {
            println!("2 - pushing {verty} onto row {i}");
        }
        log.push(verty);
        log.push('\n');
    }
    log.push(botleft);
    for i in 0..ncols {
        let mut n = maxcol[i];
        if i < ncols - 1 {
            n += sep;
        }
        for _ in 0..n {
            log.push(dash);
        }
        if vert[i] {
            if rrr[rrr.len() - 1][i + 1] != "\\ext" {
                log.push(uptee);
            } else {
                log.push(dash);
            }
            for _ in 0..sep {
                log.push(dash);
            }
        }
    }
    log.push(botright);
    log.push('\n');

    // Convert into a super-character vec of matrices.  There is one vector entry per line.
    // In each matrix, an entry is a super_character: a rust character, together with the escape
    // code characters that came before it.

    let mut mat = Vec::<Vec<Vec<char>>>::new();
    {
        let mut all = Vec::<Vec<char>>::new();
        let mut z = Vec::<char>::new();
        for c in log.chars() {
            if c != '\n' {
                z.push(c);
            } else {
                if !z.is_empty() {
                    all.push(z.clone());
                }
                z.clear();
            }
        }
        if !z.is_empty() {
            all.push(z);
        }
        for chars in all {
            mat.push(package_characters_with_escapes_char(&chars));
        }
    }

    // "Smooth" edges of hlines.

    for i in 0..mat.len() {
        for j in 0..mat[i].len() {
            if j > 0
                && mat[i][j - 1] == vec![dash]
                && mat[i][j] == vec![verty]
                && j + 1 < mat[i].len()
                && mat[i][j + 1] == vec![dash]
                && i + 1 < mat.len()
                && j < mat[i + 1].len()
                && mat[i + 1][j] != vec![verty]
            {
                mat[i][j] = vec![uptee];
            } else if j > 0
                && mat[i][j - 1] == vec![dash]
                && mat[i][j] == vec![verty]
                && j + 1 < mat[i].len()
                && mat[i][j + 1] == vec![dash]
            {
                mat[i][j] = vec![cross];
            } else if mat[i][j] == vec![verty]
                && j + 1 < mat[i].len()
                && mat[i][j + 1] == vec![dash]
            {
                mat[i][j] = vec![lefty];
            } else if j > 0
                && mat[i][j - 1] == vec![dash]
                && mat[i][j] == vec![verty]
                && (j + 1 == mat[i].len() || mat[i][j + 1] != vec![dash])
            {
                mat[i][j] = vec![righty];
            }
        }
    }

    // Output matrix.

    log.clear();
    for x in mat {
        for y in x {
            for z in y {
                log.push(z);
            }
        }
        log.push('\n');
    }

    // Finish.

    if debug_print {
        println!();
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[cfg(test)]
mod tests {

    // run this test using:
    // cargo test -p tenkit2 test_print_tabular_vbox

    use crate::print_tabular_vbox;

    // (should add some escape codes)

    #[test]
    fn test_print_tabular_vbox() {
        // test 1

        let mut rows = Vec::<Vec<String>>::new();
        let row = vec![
            "omega".to_string(),
            "superduperfineexcellent".to_string(),
            "\\ext".to_string(),
        ];
        rows.push(row);
        let row = vec![
            "woof".to_string(),
            "snarl".to_string(),
            "octopus".to_string(),
        ];
        rows.push(row);
        let row = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        rows.push(row);
        let row = vec![
            "hiccup".to_string(),
            "tomatillo".to_string(),
            "ddd".to_string(),
        ];
        rows.push(row);
        let mut log = String::new();
        let justify = &[b'r', b'|', b'l', b'l'];
        print_tabular_vbox(&mut log, &rows, 2, justify, false, false);
        let answer = "┌────────┬─────────────────────────┐\n\
                      │ omega  │  superduperfineexcellent│\n\
                      │  woof  │  snarl      octopus     │\n\
                      │     a  │  b          c           │\n\
                      │hiccup  │  tomatillo  ddd         │\n\
                      └────────┴─────────────────────────┘\n";
        if log != answer {
            println!("\ntest 1 failed");
            println!("\nyour answer:\n{log}");
            println!("correct answer:\n{answer}");
        }
        if log != answer {
            panic!();
        }

        // test 2

        let mut rows = Vec::<Vec<String>>::new();
        let row = vec!["pencil".to_string(), "pusher".to_string()];
        rows.push(row);
        let row = vec!["\\hline".to_string(), "\\hline".to_string()];
        rows.push(row);
        let row = vec!["fabulous pumpkins".to_string(), "\\ext".to_string()];
        rows.push(row);
        let mut log = String::new();
        let justify = &[b'l', b'|', b'l'];
        print_tabular_vbox(&mut log, &rows, 2, justify, false, false);
        let answer = "┌────────┬────────┐\n\
                      │pencil  │  pusher│\n\
                      ├────────┴────────┤\n\
                      │fabulous pumpkins│\n\
                      └─────────────────┘\n";
        if log != answer {
            println!("\ntest 2 failed");
            println!("\nyour answer:\n{log}");
            println!("correct answer:\n{answer}");
        }
        if log != answer {
            panic!();
        }
    }
}
