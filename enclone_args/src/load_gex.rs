// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Load gene expression and feature barcoding (antibody, antigen) data from
// Cell Ranger outputs.

use crate::load_gex_core::load_gex;
use enclone_core::defs::{EncloneControl, GexInfo};

use hdf5::Dataset;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write;
use vector_utils::{bin_position, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Get gene expression and feature barcoding counts.

pub fn get_gex_info(ctl: &mut EncloneControl) -> Result<GexInfo, String> {
    let mut gex_features = Vec::<Vec<String>>::new();
    let mut gex_barcodes = Vec::<Vec<String>>::new();
    let mut cluster = Vec::<HashMap<String, usize>>::new();
    let mut cell_type = Vec::<HashMap<String, String>>::new();
    let mut cell_type_specified = Vec::<bool>::new();
    let mut pca = Vec::<HashMap<String, Vec<f64>>>::new();
    let mut gex_mults = Vec::<f64>::new();
    let mut fb_mults = Vec::<f64>::new();
    let mut gex_cell_barcodes = Vec::<Vec<String>>::new();
    let mut have_gex = false;
    let mut have_fb = false;
    let mut h5_paths = Vec::<String>::new();
    let mut feature_metrics = Vec::<HashMap<(String, String), String>>::new();
    let mut json_metrics = Vec::<HashMap<String, f64>>::new();
    load_gex(
        ctl,
        &mut gex_features,
        &mut gex_barcodes,
        &mut cluster,
        &mut cell_type,
        &mut cell_type_specified,
        &mut pca,
        &mut gex_mults,
        &mut fb_mults,
        &mut gex_cell_barcodes,
        &mut have_gex,
        &mut have_fb,
        &mut h5_paths,
        &mut feature_metrics,
        &mut json_metrics,
    )?;
    if ctl.gen_opt.gene_scan.is_some() && !ctl.gen_opt.accept_inconsistent {
        let mut allf = gex_features.clone();
        unique_sort(&mut allf);
        if allf.len() != 1 {
            let mut msg = format!(
                "\nCurrently, SCAN requires that all datasets have identical \
                 features, and they do not.\n\
                There are {} datasets and {} feature sets after removal of \
                 duplicates.\nClassification of features sets:\n\n",
                gex_features.len(),
                allf.len()
            );
            for (f, id) in gex_features.iter().zip(ctl.origin_info.dataset_id.iter()) {
                let p = bin_position(&allf, f);
                writeln!(msg, "{id} ==> {p}").unwrap();
            }
            msg += "\n";
            return Err(msg);
        }
    }
    let mut h5_data = Vec::<Option<Dataset>>::new();
    let mut h5_indices = Vec::<Option<Dataset>>::new();
    let mut h5_indptr = Vec::<Vec<u32>>::new();

    let gex_outs = &ctl.origin_info.gex_path;
    for i in 0..ctl.origin_info.dataset_path.len() {
        if !gex_outs[i].is_empty() {
            let f = &h5_paths[i];

            let h = hdf5::File::open(f).unwrap();

            h5_data.push(Some(h.dataset("matrix/data").unwrap()));
            h5_indices.push(Some(h.dataset("matrix/indices").unwrap()));
            let indptr = h.dataset("matrix/indptr").unwrap();
            let x: Vec<u32> = indptr.as_reader().read().unwrap().to_vec();
            h5_indptr.push(x);
        } else {
            h5_data.push(None);
            h5_indices.push(None);
            h5_indptr.push(Vec::<u32>::new());
        }
    }

    fn compute_feature_id(gex_features: &[String]) -> HashMap<String, usize> {
        let mut x = HashMap::<String, usize>::new();
        for (j, f) in gex_features.iter().enumerate() {
            let ff = f.splitn(4, '\t').take(3).collect::<Vec<&str>>();
            for z in 0..2 {
                if ff[2].starts_with("Antibody") {
                    x.insert(format!("{}_ab", ff[z]), j);
                } else if ff[2].starts_with("CRISPR") {
                    x.insert(format!("{}_cr", ff[z]), j);
                } else if ff[2].starts_with("CUSTOM") {
                    x.insert(format!("{}_cu", ff[z]), j);
                } else if ff[2].starts_with("Gene") {
                    x.insert(format!("{}_g", ff[z]), j);
                } else if ff[2].starts_with("Antigen") {
                    x.insert(format!("{}_ag", ff[z]), j);
                }
            }
        }
        x
    }
    let n = gex_features.len();
    let pi = (0..n).into_par_iter();
    let mut feature_id = Vec::<HashMap<String, usize>>::new();
    pi.map(|i| compute_feature_id(&gex_features[i]))
        .collect_into_vec(&mut feature_id);
    let is_gex = gex_features
        .iter()
        .map(|g| {
            g.iter()
                .map(|f| {
                    let ff = f.split('\t').nth(2).unwrap();
                    ff.starts_with("Gene")
                })
                .collect()
        })
        .collect();

    // Answer.

    Ok(GexInfo {
        gex_features,
        gex_barcodes,
        cluster,
        cell_type,
        cell_type_specified,
        pca,
        gex_cell_barcodes,
        gex_mults,
        fb_mults,
        h5_data,
        h5_indices,
        h5_indptr,
        is_gex,
        feature_id,
        have_gex,
        have_fb,
        feature_metrics,
        json_metrics,
    })
}
