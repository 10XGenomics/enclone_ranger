// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use self::refx::RefData;
use crate::{
    barcode_fate::BarcodeFate,
    defs::{AlleleData, CloneInfo, EncloneControl, ExactClonotype, GexInfo},
};
use enclone_proto::types::DonorReferenceItem;
use qd::Double;
use std::{collections::HashMap, time::Instant};
use vdj_ann::refx;

#[derive(Clone, Debug, Default)]
pub struct MainEncloneOutput {
    pub pics: Vec<String>, // clonotype tables
    pub last_widths: Vec<u32>,
    pub svgs: Vec<String>, // SVG objects
    pub summary: String,   // summary
    pub metrics: Vec<String>,
    pub dataset_names: Vec<String>,
    pub parseable_stdouth: bool,
    pub noprint: bool,
    pub noprintx: bool,
    pub html: bool,
    pub ngroup: bool,
    pub pretty: bool,
}

#[derive(Default)]
pub struct EncloneState {
    pub inter: EncloneIntermediates,
    pub outs: MainEncloneOutput,
}

#[derive(Default)]
pub struct EncloneSetup {
    pub ctl: EncloneControl,
    pub ann: String,
    pub gex_info: GexInfo,
    pub tall: Option<Instant>,
    pub refdata: RefData,
    pub is_bcr: bool,
    pub to_ref_index: HashMap<usize, usize>,
}

#[derive(Default)]
pub struct EncloneIntermediates {
    pub setup: EncloneSetup,
    pub ex: EncloneExacts,
}

pub type BarcodeFates = HashMap<String, BarcodeFate>;

#[derive(Default, Clone)]
pub struct EncloneExacts {
    pub to_bc: HashMap<(usize, usize), Vec<String>>,
    pub exact_clonotypes: Vec<ExactClonotype>,
    pub raw_joins: Vec<Vec<usize>>,
    pub info: Vec<CloneInfo>,
    pub orbits: Vec<Vec<i32>>,
    pub vdj_cells: Vec<Vec<String>>,
    pub join_info: Vec<JoinInfo>,
    pub drefs: Vec<DonorReferenceItem>,
    pub sr: Vec<Vec<Double>>,
    pub is_bcr: bool,
    pub allele_data: AlleleData,
}

#[derive(Clone)]
pub struct JoinInfo {
    pub index1: usize,
    pub index2: usize,
    pub err: bool,
    pub log: Vec<u8>,
}
