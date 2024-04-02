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

#[derive(Default)]
pub struct EncloneSetup {
    pub ctl: EncloneControl,
    pub ann: String,
    pub gex_info: GexInfo,
    pub tall: Option<Instant>,
    pub refdata: RefData,
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
    pub allele_data: AlleleData,
}

#[derive(Clone)]
pub struct JoinInfo {
    pub index1: usize,
    pub index2: usize,
    pub err: bool,
    pub log: Vec<u8>,
}
