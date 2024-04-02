// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a special entry point for cellranger, where we know that the arguments that could
// be passed are limited.  The code here is simplified and could be further simplified.

use self::refx::{make_vdj_ref_data_core, RefData};
use crate::USING_PAGER;
use enclone::innate::species;
use enclone_args::load_gex::get_gex_info;
use enclone_args::proc_args::proc_args;
use enclone_core::defs::EncloneControl;
use enclone_core::enclone_structs::EncloneSetup;
use enclone_print::process_clonotypes::{process_clonotypes, OrbitProcessor};
use enclone_stuff::start::main_enclone_start;
use std::sync::atomic::Ordering::SeqCst;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};
use string_utils::TextUtils;
use vdj_ann::refx;

pub fn main_enclone_ranger(args: &[String]) -> Result<(), String> {
    const REQUIRED_ARGS: [&str; 8] = [
        "CELLRANGER",
        "DONOR_REF_FILE",
        "MAX_CORES",
        "NOPAGER",
        "NOPRINT",
        "PRE",
        "PROTO",
        "REF",
    ];
    const ALLOWED_ARGS: [&str; 17] = [
        "BCR",
        "META",
        "NOPRETTY",
        "PROTO_METADATA",
        "TCR",
        "TCRGD",
        "GAMMA_DELTA",
        "FATE_FILE",
        "NUMI",
        "NUMI_RATIO",
        "NGRAPH_FILTER",
        "NWEAK_CHAINS",
        "NFOURSIE_KILL",
        "NDOUBLET",
        "NSIG",
        "SPLIT_MAX_CHAINS",
        "NCROSS",
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
            panic!("Illegal argument {arg} passed to main_enclone_ranger.");
        }
    }
    for (found, arg) in found.into_iter().zip(REQUIRED_ARGS.into_iter()) {
        if !found {
            panic!("Required argument {arg} not passed to main_enclone_ranger");
        }
    }
    let setup = main_enclone_setup_ranger(args)?;
    let (exacts, fate) = main_enclone_start(&setup)?;
    let gex_readers = setup.create_gex_readers();
    process_clonotypes::<(), ()>(&setup, &exacts, &gex_readers, &fate, NoOpProc)
}

pub fn main_enclone_setup_ranger(args: &[String]) -> Result<EncloneSetup, String> {
    let tall = Instant::now();

    // Set up stuff, read args, etc.

    let mut ctl = EncloneControl::default();
    ctl.gen_opt.cellranger = true;
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
    make_vdj_ref_data_core(
        &mut refdata,
        refx2,
        &ext_refx,
        ctl.gen_opt.is_tcr(),
        ctl.gen_opt.is_bcr(),
        None,
    );

    // Determine if the species is human or mouse or unknown.

    ctl.gen_opt.species = species(&refdata).to_string();

    // Return.

    Ok(EncloneSetup {
        ctl,
        refdata,
        ann: ann.to_string(),
        gex_info,
        tall: Some(tall),
    })
}

struct NoOpProc;

impl OrbitProcessor<(), ()> for NoOpProc {}
