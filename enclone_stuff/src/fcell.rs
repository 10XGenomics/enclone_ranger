// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Filter using constraints imposed by FCELL.

use enclone_core::defs::{CloneInfo, ExactClonotype};
use enclone_core::enclone_structs::EncloneSetup;
use enclone_print::get_gex_matrix_entry::get_gex_matrix_entry;
use evalexpr::{ContextWithMutableVariables, HashMapContext};

use vector_utils::{bin_position, erase_if};

pub fn filter_by_fcell(
    setup: &EncloneSetup,
    orbits: &mut Vec<Vec<i32>>,
    info: &[CloneInfo],
    exact_clonotypes: &mut [ExactClonotype],
) -> Result<(), String> {
    let (gex_info, ctl) = (&setup.gex_info, &setup.ctl);
    if !ctl.clono_filt_opt_def.fcell.is_empty() {
        let gex_readers = setup.create_gex_readers();

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
                        if p >= 0 {
                            let z1 = gex_info.h5_indptr[li][p as usize] as usize;
                            let z2 = gex_info.h5_indptr[li][p as usize + 1] as usize; // p+1 OK?
                            (d_all[l], ind_all[l]) =
                                gex_readers[li].as_ref().unwrap().get_range(z1..z2).unwrap();
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
                    for x in &ctl.clono_filt_opt_def.fcell {
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
                                            ctl, gex_info, fid, &d_all, &ind_all, li, l, var,
                                        );
                                        val = format!("{raw_count:.2}");
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
