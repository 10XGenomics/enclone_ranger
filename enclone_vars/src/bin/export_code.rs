// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Actually export code.

// Read the vars file and export code.  This is a partial implementation.

use enclone_vars::export_code::export_code;
use io_utils::{fwrite, open_for_write_new};
use pretty_trace::PrettyTrace;

use std::io::Write;

fn main() {
    PrettyTrace::new().on();
    let outs = export_code(0);
    for out in outs {
        let mut f = open_for_write_new![&out.0];
        fwrite!(f, "{}", out.1);
    }
}
