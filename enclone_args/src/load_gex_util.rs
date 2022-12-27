// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::EncloneControl;
use io_utils::{dir_list, path_exists};
use vector_utils::VecUtils;

pub fn find_pca_file(
    _ctl: &EncloneControl,
    _outs: &str,
    analysis: &[String],
    pathlist: &mut Vec<String>,
) -> String {
    let mut pca_file = String::new();
    for x in analysis.iter() {
        pca_file = format!("{x}/pca/10_components/projection.csv");
        if path_exists(&pca_file) {
            pathlist.push(pca_file.clone());
            break;
        }
        pca_file = format!("{x}/pca/gene_expression_10_components/projection.csv");
        if path_exists(&pca_file) {
            pathlist.push(pca_file.clone());
            break;
        }
    }
    pca_file
}

pub fn find_json_metrics_file(
    ctl: &EncloneControl,
    _outs: &str,
    analysis: &[String],
    pathlist: &mut Vec<String>,
) -> String {
    let mut json_metrics_file = String::new();
    if !ctl.gen_opt.cellranger {
        for x in analysis.iter() {
            let f = format!("{x}/metrics_summary_json.json");
            if path_exists(&f) {
                json_metrics_file = f.clone();
                pathlist.push(f);
                break;
            }
        }
    }
    json_metrics_file
}

pub fn find_feature_metrics_file(
    ctl: &EncloneControl,
    _outs: &str,
    analysis: &[String],
    pathlist: &mut Vec<String>,
) -> String {
    let mut feature_metrics_file = String::new();
    if !ctl.gen_opt.cellranger {
        for x in analysis.iter() {
            let f = format!("{x}/per_feature_metrics.csv");
            if path_exists(&f) {
                feature_metrics_file = f.clone();
                pathlist.push(f);
                break;
            }
        }
    }
    feature_metrics_file
}

pub fn find_metrics_file(
    ctl: &EncloneControl,
    outs: &str,
    _analysis: &[String],
    pathlist: &mut Vec<String>,
) -> String {
    let mut metrics_file = String::new();
    if !ctl.gen_opt.cellranger {
        let summary_dir = format!("{outs}/../multi_web_summary_json/metrics_summary_csv");
        if path_exists(&summary_dir) {
            let list = dir_list(&summary_dir);
            if list.solo() {
                let path = format!("{summary_dir}/{}", list[0]);
                pathlist.push(path.clone());
                metrics_file = path;
            }
        }
    }
    metrics_file
}

pub fn find_cluster_file(
    _ctl: &EncloneControl,
    _outs: &str,
    analysis: &[String],
    pathlist: &mut Vec<String>,
) -> String {
    let mut cluster_file = String::new();
    for x in analysis.iter() {
        cluster_file = format!("{x}/clustering/graphclust/clusters.csv");
        if path_exists(&cluster_file) {
            pathlist.push(cluster_file.clone());
            break;
        }
        cluster_file = format!("{x}/clustering/gene_expression_graphclust/clusters.csv");
        if path_exists(&cluster_file) {
            pathlist.push(cluster_file.clone());
            break;
        }
    }
    cluster_file
}
