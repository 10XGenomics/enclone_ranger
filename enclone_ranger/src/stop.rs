// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::enclone_structs::{BarcodeFates, EncloneIntermediates};
use enclone_print::print_clonotypes::print_clonotypes;
use hdf5::Reader;
use rayon::prelude::*;

pub fn main_enclone_stop_ranger(
    inter: EncloneIntermediates,
    fate: Vec<BarcodeFates>,
) -> Result<(), String> {
    // Unpack inputs.
    let gex_info = &inter.setup.gex_info;
    let ctl = &inter.setup.ctl;

    // Load the GEX and FB data.  This is quite horrible: the code and computation are duplicated
    // verbatim in fcell.rs.

    let mut d_readers = Vec::<Option<Reader<'_>>>::new();
    let mut ind_readers = Vec::<Option<Reader<'_>>>::new();
    for li in 0..ctl.origin_info.n() {
        if !ctl.origin_info.gex_path[li].is_empty() {
            let x = gex_info.h5_data[li].as_ref();
            d_readers.push(Some(x.unwrap().as_reader()));
            ind_readers.push(Some(gex_info.h5_indices[li].as_ref().unwrap().as_reader()));
        } else {
            d_readers.push(None);
            ind_readers.push(None);
        }
    }
    let mut h5_data = Vec::<(usize, Vec<u32>, Vec<u32>)>::new();
    for li in 0..ctl.origin_info.n() {
        h5_data.push((li, Vec::new(), Vec::new()));
    }
    h5_data.par_iter_mut().for_each(|res| {
        let li = res.0;
        if !ctl.origin_info.gex_path[li].is_empty() && ctl.gen_opt.h5_pre {
            res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
            res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
        }
    });

    // Find and print clonotypes.  (But we don't actually print them here.)
    print_clonotypes(
        &inter.setup,
        &inter.ex,
        &d_readers,
        &ind_readers,
        &h5_data,
        &fate,
    )?;
    Ok(())
}
