// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Load gene expression and feature barcoding (antibody, antigen) data from Cell Ranger outputs.

use crate::load_gex_util::{
    find_cluster_file, find_feature_metrics_file, find_json_metrics_file, find_metrics_file,
    find_pca_file,
};
use crate::{fnx, parse_csv_pure};
use enclone_core::defs::EncloneControl;
use enclone_core::slurp::slurp_h5;
use io_utils::{dir_list, open_for_read, open_userfile_for_read, path_exists};
use itertools::Itertools;
use rayon::prelude::*;
use serde_json::Value;
use std::{
    collections::HashMap,
    convert::TryInto,
    fmt::Write,
    fs::read_to_string,
    io::{BufRead, Read},
    time::Instant,
};
use string_utils::{parse_csv, TextUtils};
use vector_utils::{unique_sort, VecUtils};

#[derive(Default)]
struct LoadResult {
    gex_features: Vec<String>,
    gex_barcodes: Vec<String>,
    gex_mult: Option<f64>,
    fb_mult: Option<f64>,
    gex_cell_barcodes: Vec<String>,
    cluster: HashMap<String, usize>,
    cell_type: HashMap<String, String>,
    pca: HashMap<String, Vec<f64>>,
    cell_type_specified: bool,
    error: String,
    h5_path: String,
    f15: Vec<String>,
    feature_metrics: HashMap<(String, String), String>,
    json_metrics: HashMap<String, f64>,
    metrics: String,
    fb_total_umis: u64,
    fb_brn: Vec<(String, u32, u32)>,
    feature_refs: String,
    fb_brnr: Vec<(String, u32, u32)>,
    fb_total_reads: u64,
    fb_bdcs: Vec<(String, u32, u32, u32)>,
}

pub fn load_gex(
    ctl: &mut EncloneControl,
    gex_features: &mut Vec<Vec<String>>,
    gex_barcodes: &mut Vec<Vec<String>>,
    fb_top_reads_barcodes: &mut Vec<Vec<String>>,
    fb_total_umis: &mut Vec<u64>,
    fb_total_reads: &mut Vec<u64>,
    fb_brn: &mut Vec<Vec<(String, u32, u32)>>,
    fb_brnr: &mut Vec<Vec<(String, u32, u32)>>,
    fb_bdcs: &mut Vec<Vec<(String, u32, u32, u32)>>,
    feature_refs: &mut Vec<String>,
    cluster: &mut Vec<HashMap<String, usize>>,
    cell_type: &mut Vec<HashMap<String, String>>,
    cell_type_specified: &mut Vec<bool>,
    pca: &mut Vec<HashMap<String, Vec<f64>>>,
    gex_mults: &mut Vec<f64>,
    fb_mults: &mut Vec<f64>,
    gex_cell_barcodes: &mut Vec<Vec<String>>,
    have_gex: &mut bool,
    have_fb: &mut bool,
    h5_paths: &mut Vec<String>,
    feature_metrics: &mut Vec<HashMap<(String, String), String>>,
    json_metrics: &mut Vec<HashMap<String, f64>>,
    metrics: &mut Vec<String>,
) -> Result<(), String> {
    let t = Instant::now();
    let mut results = Vec::<(usize, LoadResult)>::new();
    for i in 0..ctl.origin_info.gex_path.len() {
        results.push((i, LoadResult::default()));
    }
    let gex_outs = &ctl.origin_info.gex_path;
    // Here and in other places, where an error message can be printed in a parallel loop, it
    // would be better if the thread could use a global lock to prevent multiple threads from
    // issuing an error message.
    //
    // A lot of time is spent in this parallel loop.  Some things are known about this:
    // 1. When running it over a large number of datasets, the observed load average is ~2, so
    //    somehow the parallelism is not working.
    // 2. We know where the time is spent in the loop, and this is marked below.
    results.par_iter_mut().for_each(|(i, r)| {
        let pathlist = &mut r.f15;
        let i = *i;
        if !gex_outs[i].is_empty() {
            // First define the path where the GEX files should live, and make sure that the path
            // exists.

            let root = gex_outs[i].clone();
            let mut outs = root.clone();
            if root.ends_with("/outs") && path_exists(&root) {
                outs = root;
            } else if root.ends_with("/outs") {
                outs = root.before("/outs").to_string();
                if !path_exists(&outs) {
                    r.error = format!(
                        "\nThe directory\n{outs}\ndoes not exist.  Something must be amiss with \
                        the arguments to PRE and/or GEX and/or META.\n"
                    );
                    return;
                }
            }

            // Define the file paths and test for their existence.

            let mut h5_path = String::new();
            let h5p = [
                "raw_feature_bc_matrix.h5",
                "raw_gene_bc_matrices_h5.h5",
                "multi/count/raw_feature_bc_matrix.h5",
            ];
            for x in h5p.iter() {
                let p = format!("{outs}/{x}");
                if path_exists(&p) {
                    pathlist.push(p.clone());
                    h5_path = p;
                    break;
                }
            }
            if h5_path.is_empty() {
                r.error = format!(
                    "\nThe file raw_feature_bc_matrix.h5 is not in the directory\n{outs}\n\
                    and neither is the older-named version raw_gene_bc_matrices_h5.h5.  Perhaps \
                    something\nis amiss with the arguments to PRE and/or GEX and/or META.\n"
                );
                return;
            }
            r.h5_path = h5_path.clone();
            let types_file = format!("{outs}/analysis_csv/celltypes/celltypes.csv");

            // Define possible places for the analysis directory.

            let mut analysis = Vec::<String>::new();
            analysis.push(outs.to_string());
            analysis.push(format!("{outs}/analysis_csv"));
            analysis.push(format!("{outs}/analysis"));
            analysis.push(format!("{outs}/count/analysis"));
            let pso1 = format!("{outs}/per_sample_outs");
            let pso2 = format!("{outs}/../per_sample_outs");
            for pso in [pso1, pso2].iter() {
                if path_exists(pso) {
                    let samples = dir_list(pso);
                    if samples.solo() {
                        let a = format!("{pso}/{}/count/analysis", samples[0]);
                        analysis.push(a);
                        let a = format!("{pso}/{}/count/analysis_csv", samples[0]);
                        analysis.push(a);
                    }
                }
            }

            // Find files.

            let pca_file = find_pca_file(ctl, &outs, &analysis, pathlist);
            let json_metrics_file = find_json_metrics_file(ctl, &outs, &analysis, pathlist);
            let feature_metrics_file = find_feature_metrics_file(ctl, &outs, &analysis, pathlist);
            let metrics_file = find_metrics_file(ctl, &outs, &analysis, pathlist);
            let cluster_file = find_cluster_file(ctl, &outs, &analysis, pathlist);

            // Proceed.

            for f in [pca_file.clone(), cluster_file.clone()].iter() {
                if !path_exists(f) {
                    r.error = format!(
                        "\nThe file\n{f}\ndoes not exist.  \
                        Perhaps one of your directories is missing some stuff.\n\n\
                        One possibility is that you ran \"cellranger count\" using only \
                        feature barcode (antibody) data,\nand you had less then ten antibodies.  \
                        Currently if you do this, cellranger will not run the\nsecondary \
                        analyses, so you'll be missing some files.  A workaround is to add \
                        some \"fake\" antibodies\nto pad out the total number to ten.\n\n\
                        Another possibility is that this is a multi run, and the path you \
                        provided\nis to a subdirectory of the outs folder.  In that case it may \
                        work to provide the path to outs\nor (equivalently) the parent \
                        directory.\n"
                    );
                    return;
                } else {
                    pathlist.push(f.to_string());
                }
            }

            // Find metrics summary file.

            let mut csv = String::new();
            let mut csvs = Vec::<String>::new();
            csvs.push(format!("{outs}/metrics_summary.csv"));
            csvs.push(format!("{outs}/metrics_summary_csv.csv"));
            let pso = format!("{outs}/per_sample_outs");
            if path_exists(&pso) {
                let samples = dir_list(&pso);
                if samples.solo() {
                    let a = format!("{pso}/{}/metrics_summary.csv", samples[0]);
                    csvs.push(a);
                    let a = format!("{pso}/{}/metrics_summary_csv.csv", samples[0]);
                    csvs.push(a);
                }
            }
            for c in &csvs {
                if path_exists(c) {
                    csv = c.clone();
                    pathlist.push(c.to_string());
                    break;
                }
            }
            if csv.is_empty() {
                r.error = format!(
                    "\nSomething wrong with GEX or META argument:\ncan't find the file \
                        metrics_summary.csv or metrics_summary_csv.csv in the directory\n\
                        {outs}"
                );
                return;
            }

            // Read cell types.

            if path_exists(&types_file) {
                pathlist.push(types_file.clone());
                let f = open_userfile_for_read(&types_file);
                let mut count = 0;
                for line in f.lines() {
                    count += 1;
                    if count == 1 {
                        continue;
                    }
                    let s = line.unwrap();
                    let barcode = s.before(",");
                    let cell_type = s.after(",");
                    r.cell_type.insert(barcode.to_string(), cell_type.to_string());
                    r.cell_type_specified = true;
                }
            } else if ctl.gen_opt.mark_stats
                || ctl.gen_opt.mark_stats2
                || ctl.clono_filt_opt_def.marked_b
            {
                r.error = format!(
                    "\nIf you use MARK_STATS or MARK_STATS2 or MARKED_B, celltypes.csv has to \
                    exist, and this file\n{types_file}\ndoes not exist.\n"
                );
                return;
            }

            // Read json metrics file.  Note that we do not enforce the requirement of this
            // file, so it may not be present.  Also it is not present in the outs folder of CS
            // pipelines, and a customer would have to rerun with --vdrmode=disable to avoid
            // deleting the file, and then move it to outs so enclone could find it.

            if !json_metrics_file.is_empty() {
                let m = std::fs::read_to_string(&json_metrics_file).unwrap();
                let v: Value = serde_json::from_str(&m).unwrap();
                let z = v.as_object().unwrap();
                for (var, value) in z.iter() {
                    if value.as_f64().is_some() {
                        let value = value.as_f64().unwrap();
                        r.json_metrics.insert(var.to_string(), value);
                    }
                }
            }

            // Read and parse metrics file.  Rewrite as metrics class, metric name, metric value.

            if !metrics_file.is_empty() {
                let m = std::fs::read_to_string(&metrics_file).unwrap();
                let fields = parse_csv_pure(m.before("\n"));
                let (mut class, mut name, mut value) = (None, None, None);
                for field in fields {
                    if field == "Library Type" {
                        class = Some(i);
                    } else if field == "Metric Name" {
                        name = Some(i);
                    } else if field == "Metric Value" {
                        value = Some(i);
                    }
                }
                let (class, name, value) = (class.unwrap(), name.unwrap(), value.unwrap());
                let mut lines = Vec::<String>::new();
                let mut first = true;
                for line in m.lines() {
                    if first {
                        first = false;
                    } else {
                        let fields = parse_csv_pure(line);
                        lines.push(format!(
                            "{},{},{}",
                            fields[class], fields[name], fields[value]
                        ));
                    }
                }
                r.metrics = format!("{}\n", lines.iter().format("\n"));
            }

            // Read feature metrics file.  Note that we do not enforce the requirement of this
            // file, so it may not be present.

            if !feature_metrics_file.is_empty() {
                let f = open_for_read![&feature_metrics_file];
                let mut feature_pos = HashMap::<String, usize>::new();
                let mut xfields = Vec::<String>::new();
                for (count, line) in f.lines().enumerate() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    if count == 0 {
                        for (j,field) in fields.iter().enumerate() {
                            feature_pos.insert(field.to_string(), j);
                        }
                        xfields = fields.clone();
                    } else {
                        let feature_type = &fields[feature_pos["feature_type"]];
                        let mut feature;
                        for pass in 1..=2 {
                            if pass == 1 {
                                feature = fields[feature_pos["feature_name"]].clone();
                            } else {
                                feature = fields[feature_pos["feature_id"]].clone();
                            }
                            if feature_type.starts_with("Antibody") {
                                feature += "_ab";
                            } else if feature_type.starts_with("CRISPR") {
                                feature += "_cr";
                            } else if feature_type.starts_with("CUSTOM") {
                                feature += "_cu";
                            } else if feature_type.starts_with("Gene") {
                                feature += "_g";
                            } else if feature_type.starts_with("Antigen") {
                                feature += "_ag";
                            }
                            for j in 0..fields.len() {
                                if xfields[j] == "num_umis"
                                    || xfields[j] == "num_reads"
                                    || xfields[j] == "num_umis_cells"
                                    || xfields[j] == "num_reads_cells"
                                {
                                    r.feature_metrics.insert(
                                        (feature.clone(), xfields[j].clone()),
                                        fields[j].clone(),
                                    );
                                }
                            }
                        }
                    }
                }
            }

            // Read PCA file.

            let f = open_userfile_for_read(&pca_file);
            let mut count = 0;
            for line in f.lines() {
                count += 1;
                if count == 1 {
                    continue;
                }
                let s = line.unwrap();
                let barcode = s.before(",");
                let y = s.after(",").split(',').map(str::force_f64).collect();
                // This assert is turned off because in fact there are not always 10 components.
                // assert_eq!(x.len(), 10);
                r.pca.insert(barcode.to_string(), y);
            }

            // Read graph clusters, and also get the cell barcodes from that.

            let f = open_userfile_for_read(&cluster_file);
            let mut count = 0;
            for line in f.lines() {
                count += 1;
                if count == 1 {
                    continue;
                }
                let s = line.unwrap();
                let (barcode, cluster) = (s.before(","), s.after(",").force_usize());
                r.cluster.insert(barcode.to_string(), cluster);
                r.gex_cell_barcodes.push(barcode.to_string());
            }

            // Get the multipliers gene and feature barcode counts.

            let (mut gene_mult, mut fb_mult) = (None, None);
            let (mut rpc, mut fbrpc) = (None, None);
            let mut lines = Vec::<String>::new();
            {
                let f = open_userfile_for_read(&csv);
                for line in f.lines() {
                    let s = line.unwrap();
                    lines.push(s.to_string());
                }
            }
            if lines.is_empty() {
                r.error = format!("\nThe file\n{csv}\nis empty.\n");
                return;
            }
            let fields = parse_csv(&lines[0]);
            if fields.contains(&"Metric Name".to_string())
                && fields.contains(&"Metric Value".to_string())
                && fields.contains(&"Library Type".to_string())
            {
                let mut lib_field = 0;
                let mut name_field = 0;
                let mut value_field = 0;
                for (i, field) in fields.iter().enumerate() {
                    if field == "Library Type" {
                        lib_field = i;
                    } else if field == "Metric Name" {
                        name_field = i;
                    } else if field == "Metric Value" {
                        value_field = i;
                    }
                }
                for (j,line) in lines.iter().enumerate().skip(1) {
                    let fields = parse_csv(line);
                    if fields.len() < lib_field + 1
                        || fields.len() < name_field + 1
                        || fields.len() < value_field + 1
                    {
                        r.error = format!(
                            "\nSomething appears to be wrong with the file\n{}:\n\
                            line {} doesn't have enough fields.\n",
                            csv,
                            j + 1,
                        );
                        return;
                    }
                    if fields[lib_field] == "Gene Expression"
                        && fields[name_field] == "Mean reads per cell"
                    {
                        let mut rpcx = fields[value_field].to_string();
                        rpcx = rpcx.replace(',', "");
                        rpcx = rpcx.replace('\"', "");
                        if rpcx.parse::<usize>().is_err() {
                            r.error = format!(
                                "\nSomething appears to be wrong with the file\n{csv}:\n\
                                the Gene Expression Mean Reads per Cell value isn't an integer.\n"
                            );
                            return;
                        }
                        rpc = Some(rpcx.force_usize() as isize);
                    // Note that where we have "Antibody Capture"/"Antigen Capture", we could hypothetically have
                    // "CRISPR Guide Capture" or "Custom Feature".
                    } else if (fields[lib_field] == "Antibody Capture" || fields[lib_field] == "Antigen Capture")
                        && fields[name_field] == "Mean reads per cell"
                    {
                        let mut fbrpcx = fields[value_field].to_string();
                        fbrpcx = fbrpcx.replace(',', "");
                        fbrpcx = fbrpcx.replace('\"', "");
                        if fbrpcx.parse::<usize>().is_err() {
                            r.error = format!(
                                "\nSomething appears to be wrong with the file\n{csv}:\n\
                                the Antibody/Antigen Capture Mean Reads per Cell value isn't an integer.\n"
                            );
                            return;
                        }
                        fbrpc = Some(fbrpcx.force_usize() as isize);
                    }
                }
                if rpc.is_none() && fbrpc.is_none() {
                    r.error = format!(
                        "\nGene expression or feature barcode data was expected, however the \
                        CSV file\n{csv}\n\
                        does not have values for Gene Expression Mean Reads per Cell or
                        Antibody/Antigen Capture Mean Reads per Cell.\n\
                        This is puzzling.\n",
                    );
                    return;
                }
            } else {
                let (mut rpc_field, mut fbrpc_field) = (None, None);
                for (line_no,line) in lines.iter().enumerate() {
                    let s = line;
                    let fields = parse_csv(s);
                    if line_no == 0 {
                        for (i,field) in fields.iter().enumerate() {
                            if field == "Mean Reads per Cell" {
                                rpc_field = Some(i);
                            } else if field == "Antibody: Mean Reads per Cell" || field == "Antigen: Mean Reads per Cell"{
                                fbrpc_field = Some(i);
                            }
                        }
                    } else if line_no == 1 {
                        if rpc_field.is_some() && rpc_field.unwrap() >= fields.len() {
                            r.error = format!(
                                "\nSomething appears to be wrong with the file\n{csv}:\n\
                                the second line doesn't have enough fields.\n"
                            );
                            return;
                        } else if rpc_field.is_some() {
                            let mut rpcx = fields[rpc_field.unwrap()].to_string();
                            rpcx = rpcx.replace(',', "");
                            rpcx = rpcx.replace('\"', "");
                            if rpcx.parse::<usize>().is_err() {
                                r.error = format!(
                                    "\nSomething appears to be wrong with the file\n{csv}:\n\
                                    the Mean Reads per Cell field isn't an integer.\n"
                                );
                                return;
                            }
                            rpc = Some(rpcx.force_usize() as isize);
                        }
                        if fbrpc_field.is_some() && fbrpc_field.unwrap() >= fields.len() {
                            r.error = format!(
                                "\nSomething appears to be wrong with the file\n{csv}:\n\
                                the second line doesn't have enough fields.\n"
                            );
                            return;
                        } else if fbrpc_field.is_some() {
                            let mut fbrpcx = fields[fbrpc_field.unwrap()].to_string();
                            fbrpcx = fbrpcx.replace(',', "");
                            fbrpcx = fbrpcx.replace('\"', "");
                            if fbrpcx.parse::<usize>().is_err() {
                                r.error = format!(
                                    "\nSomething appears to be wrong with the file\n{csv}:\n\
                                    the Antibody/Antigen: Mean Reads per Cell field isn't an integer.\n"
                                );
                                return;
                            }
                            fbrpc = Some(fbrpcx.force_usize() as isize);
                        }
                    }
                }
                if rpc.is_none() && fbrpc.is_none() {
                    r.error = format!(
                        "\nGene expression or feature barcode data was expected, however the \
                        CSV file\n{csv}\n\
                        does not have a field \"Mean Reads per Cell\" or \
                        \"Antibody: Mean Reads per Cell\".\n\
                        This is puzzling, and might be because a file within the Cell Ranger outs \
                        directory has been moved\n\
                        from its original location.\n",
                    );
                    return;
                }
            }
            if let Some(rpc) = rpc {
                const RPC_EXPECTED: f64 = 20_000.0;
                gene_mult = Some(RPC_EXPECTED / rpc as f64);
            }
            if let Some(fbrpc) = fbrpc {
                const FB_RPC_EXPECTED: f64 = 5_000.0;
                fb_mult = Some(FB_RPC_EXPECTED / fbrpc as f64);
            }
            r.gex_mult = gene_mult;
            r.fb_mult = fb_mult;

            // Read the total UMIs.

            let top_file = fnx(&outs, "feature_barcode_matrix_top.total");
            if path_exists(&top_file) {
                pathlist.push(top_file.clone());
                let mut f = open_for_read![&top_file];
                let mut bytes = Vec::<u8>::new();
                f.read_to_end(&mut bytes).unwrap();
                r.fb_total_umis = u64::from_ne_bytes(bytes.try_into().unwrap());
            }

            // Read the total reads.

            let top_file = fnx(&outs, "feature_barcode_matrix_top.total_reads");
            if path_exists(&top_file) {
                pathlist.push(top_file.clone());
                let mut f = open_for_read![&top_file];
                let mut bytes = Vec::<u8>::new();
                f.read_to_end(&mut bytes).unwrap();
                r.fb_total_reads = u64::from_ne_bytes(bytes.try_into().unwrap());
            }

            // Read the barcode-ref-nonref UMI count file.

            let brn_file = fnx(&outs, "feature_barcode_matrix_top.brn");
            if path_exists(&brn_file) {
                pathlist.push(brn_file.clone());
                let f = open_for_read![&brn_file];
                for line in f.lines() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    r.fb_brn.push((
                        fields[0].to_string(),
                        fields[1].parse::<u32>().unwrap(),
                        fields[2].parse::<u32>().unwrap(),
                    ));
                }
            }

            // Read the barcode-ref-nonref read count file.

            let brnr_file = fnx(&outs, "feature_barcode_matrix_top.brnr");
            if path_exists(&brnr_file) {
                pathlist.push(brnr_file.clone());
                let f = open_for_read![&brnr_file];
                for line in f.lines() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    r.fb_brnr.push((
                        fields[0].to_string(),
                        fields[1].parse::<u32>().unwrap(),
                        fields[2].parse::<u32>().unwrap(),
                    ));
                }
            }

            // Read the bdcs read count file.

            let bdcs_file = fnx(&outs, "feature_barcode_matrix_top.bdcs");
            if path_exists(&bdcs_file) {
                pathlist.push(bdcs_file.clone());
                let f = open_for_read![&bdcs_file];
                for line in f.lines() {
                    let s = line.unwrap();
                    let fields = parse_csv(&s);
                    r.fb_bdcs.push((
                        fields[0].to_string(),
                        fields[1].parse::<u32>().unwrap(),
                        fields[2].parse::<u32>().unwrap(),
                        fields[3].parse::<u32>().unwrap(),
                    ));
                }
            }

            // Read the feature reference file.

            let fref_file = fnx(&outs, "feature_reference.csv");
            if path_exists(&fref_file) {
                pathlist.push(fref_file.clone());
                r.feature_refs = read_to_string(&fref_file).unwrap();
            }

            // Read the feature barcode matrix file.
            if let Err(err) = slurp_h5(
                &h5_path,
                &mut r.gex_barcodes,
                &mut r.gex_features,
            ) {
                r.error = err;
                return;
            }
        }
        unique_sort(&mut r.gex_cell_barcodes);
    });
    for (_, r) in &results {
        ctl.pathlist.extend(r.f15.iter().cloned());
    }
    ctl.perf_stats(&t, "in load_gex main loop");

    // Test for error.

    let t = Instant::now();
    for (_, r) in &results {
        if !r.error.is_empty() {
            return Err(r.error.clone());
        }
    }

    // Set have_gex and have_fb.

    for (_, r) in &results {
        if r.gex_mult.is_some() {
            *have_gex = true;
        }
        if r.fb_mult.is_some() {
            *have_fb = true;
        }
    }
    h5_paths.extend(results.iter().map(|(_, r)| r.h5_path.clone()));

    // Add some metrics.

    let extras = [
        (
            "ANTIBODY_G_perfect_homopolymer_frac",
            "Antibody Capture,G Homopolymer Frac",
        ),
        (
            "GRCh38_raw_rpc_20000_subsampled_filtered_bcs_median_unique_genes_detected",
            "Gene Expression,GRCh38 Median genes per cell (20k raw reads per cell)",
        ),
        (
            "GRCh38_raw_rpc_20000_subsampled_filtered_bcs_median_counts",
            "Gene Expression,GRCh38 Median UMI counts per cell (20k raw reads per cell)",
        ),
    ];
    for x in extras.iter() {
        let metric_name = x.0.to_string();
        let metric_display_name = x.1.to_string();
        let mut have = false;
        for (_, result) in &results {
            if result.json_metrics.contains_key(&metric_name) {
                have = true;
            }
        }
        if have {
            for (_, result) in results.iter_mut() {
                let mut value = String::new();
                if result.json_metrics.contains_key(&metric_name) {
                    value = format!("{:.3}", result.json_metrics[&metric_name]);
                }
                writeln!(result.metrics, "{metric_display_name},{value}").unwrap();
            }
        }
    }

    for (_, r) in results.into_iter() {
        gex_features.push(r.gex_features);
        gex_barcodes.push(r.gex_barcodes);
        gex_mults.push(r.gex_mult.unwrap_or(1.0));
        fb_mults.push(r.fb_mult.unwrap_or(1.0));
        gex_cell_barcodes.push(r.gex_cell_barcodes);
        cluster.push(r.cluster);
        cell_type.push(r.cell_type);
        pca.push(r.pca);
        cell_type_specified.push(r.cell_type_specified);
        feature_metrics.push(r.feature_metrics);
        json_metrics.push(r.json_metrics);
        metrics.push(r.metrics);
        fb_total_umis.push(r.fb_total_umis);
        fb_brn.push(r.fb_brn);
        feature_refs.push(r.feature_refs);
        fb_brnr.push(r.fb_brnr);
        fb_total_reads.push(r.fb_total_reads);
        fb_bdcs.push(r.fb_bdcs);
    }

    // Done.

    ctl.perf_stats(&t, "in load_gex tail");
    Ok(())
}
