// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Filter using constraints imposed by FCELL.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype, GexInfo};
use enclone_print::print_utils4::get_gex_matrix_entry;
use evalexpr::{ContextWithMutableVariables, HashMapContext};
use hdf5::Reader;
use io_utils::{dir_list, path_exists};
use ndarray::s;
use rayon::prelude::*;
use std::env;
use std::thread;
use std::time;
use std::time::Instant;
use vector_utils::{bin_position, erase_if};

pub fn filter_by_fcell(
    ctl: &EncloneControl,
    orbits: &mut Vec<Vec<i32>>,
    info: &[CloneInfo],
    exact_clonotypes: &mut [ExactClonotype],
    gex_info: &GexInfo,
) -> Result<(), String> {
    if !ctl.clono_filt_opt_def.fcell.is_empty() {
        // Load the GEX and FB data.  This is quite horrible: the code and computation are
        // duplicated verbatim in stop.rs.

        let tdi = Instant::now();
        let mut d_readers = Vec::<Option<Reader>>::new();
        let mut ind_readers = Vec::<Option<Reader>>::new();
        for li in 0..ctl.origin_info.n() {
            if !ctl.origin_info.gex_path[li].is_empty() && !gex_info.gex_matrices[li].initialized()
            {
                let x = gex_info.h5_data[li].as_ref();
                if x.is_none() {
                    // THIS FAILS SPORADICALLY, OBSERVED MULTIPLE TIMES,
                    // CAUSING PUSH TO D_READERS BELOW TO FAIL.
                    eprintln!("\nWeird, gex_info.h5_data[li].as_ref() is None.");
                    eprintln!("Path = {}.", ctl.origin_info.gex_path[li]);
                    let current = env::current_dir().unwrap();
                    println!(
                        "The current working directory is {}",
                        current.canonicalize().unwrap().display()
                    );
                    if path_exists(&ctl.origin_info.gex_path[li]) {
                        eprintln!(
                            "The directory that is supposed to contain \
                            raw_feature_bc_matrix.h5 exists."
                        );
                        let list = dir_list(&ctl.origin_info.gex_path[li]);
                        eprintln!(
                            "This directory is {} and its contents are:",
                            ctl.origin_info.gex_path[li]
                        );
                        for (i, li) in list.into_iter().enumerate() {
                            eprintln!("{}.  {}", i + 1, li);
                        }
                        let h5_path =
                            format!("{}/raw_feature_bc_matrix.h5", ctl.origin_info.gex_path[li]);
                        eprintln!("H5 path = {}.", h5_path);
                        if !path_exists(&h5_path) {
                            let mut msg = format!("H5 path {} does not exist.\n", h5_path);
                            msg += "Retrying a few times to see if it appears.\n";
                            for _ in 0..5 {
                                msg += "Sleeping for 0.1 seconds.";
                                thread::sleep(time::Duration::from_millis(100));
                                if !path_exists(&h5_path) {
                                    msg += "Now h5 path does not exist.\n";
                                } else {
                                    msg += "Now h5 path exists.\n";
                                    break;
                                }
                            }
                            msg += "Aborting.\n";
                            return Err(msg);
                        } else {
                            println!("h5 path exists.");
                        }
                    } else {
                        println!("Path exists.");
                    }
                    println!();
                }
                d_readers.push(Some(x.unwrap().as_reader()));
                ind_readers.push(Some(gex_info.h5_indices[li].as_ref().unwrap().as_reader()));
            } else {
                d_readers.push(None);
                ind_readers.push(None);
            }
        }
        let mut h5_data = Vec::<(usize, Vec<u32>, Vec<u32>)>::new();
        for li in 0..ctl.origin_info.n() {
            h5_data.push((li, Vec::new(), Vec::new()));
        }
        h5_data.par_iter_mut().for_each(|res| {
            let li = res.0;
            if !ctl.origin_info.gex_path[li].is_empty()
                && !gex_info.gex_matrices[li].initialized()
                && ctl.gen_opt.h5_pre
            {
                res.1 = d_readers[li].as_ref().unwrap().read_raw().unwrap();
                res.2 = ind_readers[li].as_ref().unwrap().read_raw().unwrap();
            }
        });
        ctl.perf_stats(&tdi, "setting up readers, zero");

        // Proceed.

        let mut orbits2 = Vec::<Vec<i32>>::new();
        for o in orbits.iter() {
            let mut o = o.clone();
            let mut to_deletex = vec![false; o.len()];
            for j in 0..o.len() {
                let x: &CloneInfo = &info[o[j] as usize];
                let ex = &mut exact_clonotypes[x.clonotype_index];
                let mut to_delete = vec![false; ex.ncells()];
                let mut d_all = vec![Vec::<u32>::new(); ex.clones.len()];
                let mut ind_all = vec![Vec::<u32>::new(); ex.clones.len()];
                for l in 0..ex.clones.len() {
                    let li = ex.clones[l][0].dataset_index;
                    let bc = ex.clones[l][0].barcode.clone();
                    if !gex_info.gex_barcodes.is_empty() {
                        let p = bin_position(&gex_info.gex_barcodes[li], &bc);
                        if p >= 0 && !gex_info.gex_matrices[li].initialized() {
                            let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                            let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // p+1 OK?
                            if ctl.gen_opt.h5_pre {
                                d_all[l] = h5_data[li].1[z1..z2].to_vec();
                                ind_all[l] = h5_data[li].2[z1..z2].to_vec();
                            } else {
                                d_all[l] = d_readers[li]
                                    .as_ref()
                                    .unwrap()
                                    .read_slice(s![z1..z2])
                                    .unwrap()
                                    .to_vec();
                                ind_all[l] = ind_readers[li]
                                    .as_ref()
                                    .unwrap()
                                    .read_slice(s![z1..z2])
                                    .unwrap()
                                    .to_vec();
                            }
                        }
                    }
                }
                for (l, (clone, d)) in ex
                    .clones
                    .iter()
                    .take(ex.ncells())
                    .zip(to_delete.iter_mut())
                    .enumerate()
                {
                    let li = clone[0].dataset_index;
                    let bc = &clone[0].barcode;
                    let mut keep = true;
                    for x in ctl.clono_filt_opt_def.fcell.iter() {
                        let alt = &ctl.origin_info.alt_bc_fields[li];
                        let vars = x.iter_variable_identifiers().collect::<Vec<&str>>();
                        let mut vals = Vec::<String>::new();
                        for &var in &vars {
                            let mut val = String::new();
                            let mut found = false;
                            'uloop: for au in alt {
                                if au.0 == var {
                                    if let Some(v) = au.1.get(bc) {
                                        val = v.clone();
                                        found = true;
                                        break 'uloop;
                                    }
                                }
                            }
                            if !found {
                                if let Some(&fid) = gex_info.feature_id[li].get(&var.to_string()) {
                                    let p = bin_position(&gex_info.gex_barcodes[li], bc);
                                    if p >= 0 {
                                        let raw_count = get_gex_matrix_entry(
                                            ctl, gex_info, fid, &d_all, &ind_all, li, l,
                                            p as usize, var,
                                        );
                                        val = format!("{:.2}", raw_count);
                                    }
                                }
                            }
                            vals.push(val);
                        }
                        let mut c = HashMapContext::new();
                        for (&var, val) in vars.iter().zip(vals.iter()) {
                            if let Ok(val) = val.parse::<i64>() {
                                c.set_value(var.into(), evalexpr::Value::from(val)).unwrap();
                            } else if let Ok(val) = val.parse::<f64>() {
                                c.set_value(var.into(), evalexpr::Value::from(val)).unwrap();
                            } else {
                                c.set_value(var.into(), val.clone().into()).unwrap();
                            }
                        }
                        let res = x.eval_with_context(&c);
                        let ok = res == Ok(evalexpr::Value::from(true));
                        if !ok {
                            keep = false;
                        }
                    }
                    if !keep {
                        *d = true;
                    }
                }
                erase_if(&mut ex.clones, &to_delete);
                if ex.ncells() == 0 {
                    to_deletex[j] = true;
                }
            }
            erase_if(&mut o, &to_deletex);
            if !o.is_empty() {
                orbits2.push(o.clone());
            }
        }
        *orbits = orbits2;
    }
    Ok(())
}
