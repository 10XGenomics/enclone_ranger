// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

use io_utils::{dir_list, path_exists};
use vector_utils::VecUtils;

pub fn find_pca_file(analysis: &[String]) -> String {
    let mut pca_file = String::new();
    for x in analysis {
        pca_file = format!("{x}/pca/10_components/projection.csv");
        if path_exists(&pca_file) {
            break;
        }
        pca_file = format!("{x}/pca/gene_expression_10_components/projection.csv");
        if path_exists(&pca_file) {
            break;
        }
    }
    pca_file
}

pub fn find_json_metrics_file(analysis: &[String]) -> String {
    let mut json_metrics_file = String::new();
    for x in analysis {
        let f = format!("{x}/metrics_summary_json.json");
        if path_exists(&f) {
            json_metrics_file = f.clone();
            break;
        }
    }

    json_metrics_file
}

pub fn find_feature_metrics_file(analysis: &[String]) -> String {
    let mut feature_metrics_file = String::new();
    for x in analysis {
        let f = format!("{x}/per_feature_metrics.csv");
        if path_exists(&f) {
            feature_metrics_file = f.clone();
            break;
        }
    }

    feature_metrics_file
}

pub fn find_metrics_file(outs: &str) -> String {
    let mut metrics_file = String::new();
    let summary_dir = format!("{outs}/../multi_web_summary_json/metrics_summary_csv");
    if path_exists(&summary_dir) {
        let list = dir_list(&summary_dir);
        if list.solo() {
            let path = format!("{summary_dir}/{}", list[0]);
            metrics_file = path;
        }
    }

    metrics_file
}

pub fn find_cluster_file(analysis: &[String]) -> String {
    let mut cluster_file = String::new();
    for x in analysis {
        cluster_file = format!("{x}/clustering/graphclust/clusters.csv");
        if path_exists(&cluster_file) {
            break;
        }
        cluster_file = format!("{x}/clustering/gene_expression_graphclust/clusters.csv");
        if path_exists(&cluster_file) {
            break;
        }
    }
    cluster_file
}
