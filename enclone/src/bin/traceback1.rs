// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Check that an out-of-range reference within a rayon parallel loop yields a correct traceback.
// Make sure that line 17 stays as line 17.  Otherwise change the reference to traceback1.rs:17.
//
// This was originally engineered without PrettyTrace, but the problem with this was that if the
// test failed, you get a godawful mess that is impossible to distangle.

use pretty_trace::PrettyTrace;
use rayon::prelude::*;

fn main() {
    let z = vec![0; 100];
    let mut x = vec![0; 100];
    x.par_iter_mut().for_each(|r| {
        let _ = z[100 + *r]; // line 17!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    });
}

#[test]
fn test_traceback1() {
    use assert_cmd;
    use assert_cmd::prelude::*;
    use std::fmt::Write;
    use std::{env, process::Command};
    let mut cmd = Command::cargo_bin("traceback1").expect(
        "\nAttempt to run traceback1 failed.  The most likely explanation for this is that\n\
        somehow you did not run \"cargo b\".  Please try that now, and be sure you are doing\n\
        it from the top-level enclone directory.\n",
    );
    let cmd = cmd
        .output()
        .unwrap_or_else(|_| panic!("{}", "very strange, failed to execute test_traceback1"));
    let morsel = "traceback1.rs:17";
    let err = std::str::from_utf8(&cmd.stderr).unwrap();
    if !err.contains(morsel) {
        let mut head = String::new();
        let lines = err.split('\n').collect::<Vec<&str>>();
        const MAX_LINES: usize = 60;
        for &line in lines.iter().take(MAX_LINES) {
            writeln!(head, "{line}").unwrap();
        }
        eprint!(
            "\n▓▓▓ test_traceback1 failed because did not find {} as expected;\n\n\
             this was using enclone version {} \n\n\
             ▓▓▓ traceback begins with\n{}",
            morsel,
            env!("CARGO_PKG_VERSION"),
            head,
        );
        std::process::exit(1);
    }
}
