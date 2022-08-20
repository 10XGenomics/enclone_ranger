// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a special entry point for cellranger, where we know that the arguments that could
// be passed are limited.  The code here is simplified and could be further simplified.

use self::refx::{make_vdj_ref_data_core, RefData};
use crate::stop::main_enclone_stop_ranger;
use crate::USING_PAGER;
use enclone::innate::species;
use enclone_args::load_gex::get_gex_info;
use enclone_args::proc_args::proc_args;
use enclone_core::defs::EncloneControl;
use enclone_core::enclone_structs::EncloneSetup;
use enclone_stuff::start::main_enclone_start;
use std::sync::atomic::Ordering::SeqCst;
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};
use string_utils::TextUtils;
use vdj_ann::refx;

pub fn main_enclone_ranger(args: &[String]) -> Result<(), String> {
    const REQUIRED_ARGS: [&str; 9] = [
        "CELLRANGER",
        "DONOR_REF_FILE",
        "FORCE_EXTERNAL",
        "MAX_CORES",
        "NOPAGER",
        "NOPRINT",
        "PRE",
        "PROTO",
        "REF",
    ];
    const ALLOWED_ARGS: [&str; 8] = [
        "BCR",
        "META",
        "NOPRETTY",
        "PROTO_METADATA",
        "TCR",
        "TCRGD",
        "GAMMA_DELTA",
        "FATE_FILE",
    ];
    let mut found = [false; REQUIRED_ARGS.len()];
    for arg in args.iter().skip(1) {
        let mut arg = arg.as_str();
        if arg.contains('=') {
            arg = arg.before("=");
        }
        let mut ok = false;
        for (f, &x) in found.iter_mut().zip(REQUIRED_ARGS.iter()) {
            if arg == x {
                ok = true;
                *f = true;
            }
        }
        ok = ok || ALLOWED_ARGS.contains(&arg);
        if !ok {
            panic!("Illegal argument {} passed to main_enclone_ranger.", arg);
        }
    }
    for (found, arg) in found.into_iter().zip(REQUIRED_ARGS.into_iter()) {
        if !found {
            panic!(
                "Required argument {} not passed to main_enclone_ranger",
                arg
            );
        }
    }
    let setup = main_enclone_setup_ranger(args)?;
    let inter = main_enclone_start(setup)?;
    main_enclone_stop_ranger(inter)
}

pub fn main_enclone_setup_ranger(args: &[String]) -> Result<EncloneSetup, String> {
    let tall = Instant::now();

    // Set up stuff, read args, etc.

    let mut ctl = EncloneControl::default();
    ctl.gen_opt.cellranger = true;
    ctl.gen_opt.internal_run = false;
    for arg in args.iter().skip(1) {
        if arg.starts_with("PRE=") {
            ctl.gen_opt.pre.clear();
            ctl.gen_opt
                .pre
                .extend(arg.after("PRE=").split(',').map(str::to_string));
        }
    }
    ctl.start_time = Some(tall);
    ctl.gen_opt.cpu_all_start = 0;
    ctl.gen_opt.cpu_this_start = 0;
    ctl.gen_opt.nopager = true;
    ctl.pretty = true;
    ctl.gen_opt.h5 = true;
    USING_PAGER.store(false, SeqCst);
    proc_args(&mut ctl, args)?;

    // Get gene expression and feature barcode counts.

    let gex_info = get_gex_info(&mut ctl)?;

    // Determine the reference sequence that is to be used.

    let mut refx = String::new();
    let ann = "contig_annotations.json";
    let fx = File::open(&ctl.gen_opt.refname);
    let f = BufReader::new(fx.unwrap());
    for line in f.lines() {
        let s = line.unwrap();
        refx += &s;
        refx += "\n";
    }

    // Build reference data.

    let refx2 = &refx;
    let mut refdata = RefData::new();
    let ext_refx = String::new();
    let (mut is_tcr, mut is_bcr) = (true, true);
    if ctl.gen_opt.tcr {
        is_bcr = false;
    }
    if ctl.gen_opt.bcr {
        is_tcr = false;
    }
    make_vdj_ref_data_core(&mut refdata, refx2, &ext_refx, is_tcr, is_bcr, None);
    let to_ref_index: HashMap<usize, usize> = refdata
        .id
        .iter()
        .take(refdata.refs.len())
        .enumerate()
        .map(|(i, &id)| (id as usize, i))
        .collect();

    // Determine if the species is human or mouse or unknown.

    ctl.gen_opt.species = species(&refdata).to_string();

    // Return.

    Ok(EncloneSetup {
        ctl,
        refdata,
        ann: ann.to_string(),
        gex_info,
        tall: Some(tall),
        is_bcr,
        to_ref_index,
    })
}
