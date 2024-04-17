// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This is a special entry point for cellranger, where we know that the arguments that could
// be passed are limited.  The code here is simplified and could be further simplified.

use self::refx::{make_vdj_ref_data_core, RefData};
use crate::USING_PAGER;
use anyhow::anyhow;
use enclone::innate::species;
use enclone_args::load_gex::get_gex_info;
use enclone_args::proc_args::proc_args;
use enclone_core::defs::{CellrangerOpt, EncloneControl};
use enclone_core::enclone_structs::EncloneSetup;
use enclone_process::process_clonotypes::{process_clonotypes, OrbitProcessor};
use enclone_stuff::start::main_enclone_start;
use std::sync::atomic::Ordering::SeqCst;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    time::Instant,
};
use string_utils::TextUtils;
use vdj_ann::refx;

pub fn main_enclone_ranger(args: Vec<String>) -> anyhow::Result<()> {
    const REQUIRED_ARGS: [&str; 7] = [
        "CELLRANGER",     // done
        "DONOR_REF_FILE", // done
        "MAX_CORES", // FIXME: move this behavior into enclone and set thread count in CR when calling
        "NOPRINT",   // now unused in enclone_ranger
        "PRE",       // done
        "PROTO",     // done
        "REF",       // done
    ];
    const ALLOWED_ARGS: [&str; 16] = [
        "BCR",
        "META",
        "PROTO_METADATA", // done
        "TCR",
        "TCRGD",
        "GAMMA_DELTA",   // done
        "FATE_FILE",     // done
        "NUMI",          // done
        "NUMI_RATIO",    // done
        "NGRAPH_FILTER", // done
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
    let (exacts, fate) = main_enclone_start(&setup).map_err(|e| anyhow!(e))?;
    let gex_readers = setup.create_gex_readers();
    process_clonotypes::<(), ()>(&setup, &exacts, &gex_readers, &fate, NoOpProc)
        .map_err(|e| anyhow!(e))?;
    Ok(())
}

pub fn main_enclone_setup_ranger(args: Vec<String>) -> anyhow::Result<EncloneSetup> {
    let tall = Instant::now();
    // Set up stuff, read args, etc.
    let (cr_opt, args) = CellrangerOpt::from_args(args)?;

    let mut ctl = EncloneControl {
        cr_opt,
        ..Default::default()
    };

    ctl.start_time = Some(tall);
    ctl.gen_opt.cpu_all_start = 0;
    ctl.gen_opt.cpu_this_start = 0;
    ctl.pretty = true;
    USING_PAGER.store(false, SeqCst);
    proc_args(&mut ctl, &args).map_err(|e| anyhow!(e))?;

    // Get gene expression and feature barcode counts.

    let gex_info = get_gex_info(&mut ctl).map_err(|e| anyhow!(e))?;

    // Determine the reference sequence that is to be used.

    let mut refx = String::new();
    let ann = "contig_annotations.json";
    let fx = File::open(&ctl.cr_opt.refname);
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
