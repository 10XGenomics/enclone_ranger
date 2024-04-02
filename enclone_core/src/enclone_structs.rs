// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use self::refx::RefData;
use crate::{
    barcode_fate::BarcodeFate,
    defs::{AlleleData, CloneInfo, EncloneControl, ExactClonotype, GexInfo},
    h5_lazy::LazyH5Array,
};
use anyhow::Result;
use enclone_proto::types::DonorReferenceItem;
use qd::Double;
use rayon::prelude::*;
use std::{collections::HashMap, ops::Range, time::Instant};
use vdj_ann::refx;

#[derive(Default)]
pub struct EncloneSetup {
    pub ctl: EncloneControl,
    pub ann: String,
    pub gex_info: GexInfo,
    pub tall: Option<Instant>,
    pub refdata: RefData,
}

impl EncloneSetup {
    pub fn create_gex_readers(&self) -> Vec<Option<GexReaders<'_>>> {
        (0..self.ctl.origin_info.n())
            .into_par_iter()
            .map(|li| {
                if !self.ctl.origin_info.gex_path[li].is_empty() {
                    Some(GexReaders {
                        d: LazyH5Array::new(
                            self.gex_info.h5_data[li].as_ref().unwrap(),
                            self.ctl.gen_opt.h5_pre,
                        )
                        .unwrap(),
                        ind: LazyH5Array::new(
                            self.gex_info.h5_indices[li].as_ref().unwrap(),
                            self.ctl.gen_opt.h5_pre,
                        )
                        .unwrap(),
                    })
                } else {
                    None
                }
            })
            .collect()
    }
}

pub struct GexReaders<'a> {
    d: LazyH5Array<'a>,
    ind: LazyH5Array<'a>,
}

impl<'a> GexReaders<'a> {
    pub fn get_range(&self, range: Range<usize>) -> Result<(Vec<u32>, Vec<u32>)> {
        Ok((self.d.get_range(range.clone())?, self.ind.get_range(range)?))
    }
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
