// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::enclone_structs::{BarcodeFates, EncloneExacts, EncloneSetup};
use enclone_print::print_clonotypes::print_clonotypes;

pub fn main_enclone_stop_ranger(
    setup: &EncloneSetup,
    exacts: &EncloneExacts,
    fate: Vec<BarcodeFates>,
) -> Result<(), String> {
    let gex_readers = setup.create_gex_readers();

    // Find and print clonotypes.  (But we don't actually print them here.)
    print_clonotypes(setup, exacts, &gex_readers, &fate)?;
    Ok(())
}
