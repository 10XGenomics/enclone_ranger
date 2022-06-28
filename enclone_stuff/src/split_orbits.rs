// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

use enclone_core::defs::{CloneInfo, EncloneControl, ExactClonotype};
use enclone_print::define_mat::define_mat;
use enclone_proto::types::DonorReferenceItem;
use equiv::EquivRel;
use qd::Double;
use std::collections::HashMap;
use vdj_ann::refx::RefData;
use vector_utils::{bin_position, next_diff12_3, unique_sort, VecUtils};

// Check for disjoint orbits.

pub fn split_orbits(
    orbits: &mut Vec<Vec<i32>>,
    is_bcr: bool,
    to_bc: &HashMap<(usize, usize), Vec<String>>,
    sr: &[Vec<Double>],
    ctl: &EncloneControl,
    exact_clonotypes: &[ExactClonotype],
    info: &[CloneInfo],
    raw_joins: &[Vec<usize>],
    refdata: &RefData,
    dref: &[DonorReferenceItem],
) {
    let mut orbits2 = Vec::<Vec<i32>>::new();
    for o in orbits.iter() {
        let mut od = Vec::<(Vec<usize>, usize, i32)>::new();
        for id in o.iter() {
            let x: &CloneInfo = &info[*id as usize];
            od.push((x.origin.clone(), x.clonotype_id, *id));
        }
        od.sort();
        let mut exacts = Vec::<usize>::new();
        let mut j = 0;
        while j < od.len() {
            let k = next_diff12_3(&od, j as i32) as usize;
            exacts.push(od[j].1);
            j = k;
        }
        let mat = define_mat(
            is_bcr,
            to_bc,
            sr,
            ctl,
            exact_clonotypes,
            &exacts,
            &od,
            info,
            raw_joins,
            refdata,
            dref,
        );
        let cols = mat.len();

        // Define map of indices into exacts.

        let nexacts = exacts.len();
        let mut to_exacts = HashMap::<usize, usize>::with_capacity(nexacts);
        for (u, &e) in exacts.iter().enumerate() {
            to_exacts.insert(e, u);
        }

        // Get the info indices corresponding to this clonotype.

        let mut infos = Vec::<usize>::new();
        for &oi in o {
            infos.push(oi as usize);
        }

        // Define map of exacts to infos.

        let mut to_infos = vec![Vec::<usize>::new(); nexacts];
        for (i, &infoi) in infos.iter().enumerate() {
            let u = to_exacts[&info[infoi].clonotype_index];
            to_infos[u].push(i);
        }

        // Determine which columns are "left", meaning IGH or TRB.

        let mut left = vec![false; cols];
        for m in 0..cols {
            for u in 0..mat[0].len() {
                if mat[m][u].is_some() {
                    let c = mat[m][u].unwrap();
                    let ex = &exact_clonotypes[exacts[u]];
                    if ex.share[c].left {
                        left[m] = true;
                    }
                    break;
                }
            }
        }

        // Determine which pairs of configurations share both chain types, and if so, call
        // them joined.

        let mut matu = Vec::<Vec<Option<usize>>>::with_capacity(nexacts);
        for u in 0..nexacts {
            let mut m = Vec::<Option<usize>>::with_capacity(cols);
            for mm in mat.iter().take(cols) {
                m.push(mm[u]);
            }
            matu.push(m);
        }
        unique_sort(&mut matu);
        let mut eqm = vec![vec![false; matu.len()]; matu.len()];
        for (mj1, eqm) in matu.iter().zip(eqm.iter_mut()) {
            for (mj2, eqm) in matu.iter().zip(eqm.iter_mut()) {
                let (mut l, mut r) = (false, false);
                for ((&mm1, &mm2), &ll) in mj1.iter().zip(mj2.iter()).zip(left.iter()).take(cols) {
                    if mm1.is_some() && mm2.is_some() {
                        if ll {
                            l = true;
                        } else {
                            r = true;
                        }
                    }
                }
                if l && r {
                    *eqm = true;
                }
            }
        }

        // Propagate this to an equivalence relation on the orbit elements.

        let mut eqx = EquivRel::new(o.len() as i32);
        let mut lists = vec![Vec::<usize>::new(); matu.len()];
        for u in 0..nexacts {
            let mut m = Vec::<Option<usize>>::with_capacity(cols);
            for mat in mat.iter().take(cols) {
                m.push(mat[u]);
            }
            lists[bin_position(&matu, &m) as usize].push(u);
        }
        for (l1, eqm) in lists.iter().zip(eqm.into_iter()) {
            for (l2, eqm) in lists.iter().zip(eqm.into_iter()) {
                if eqm {
                    let u1 = l1[0];
                    for &u2 in l2.iter() {
                        for &i1 in to_infos[u1].iter() {
                            for &i2 in to_infos[u2].iter() {
                                eqx.join(i1 as i32, i2 as i32);
                            }
                        }
                    }
                }
            }
        }
        let mut reps = Vec::<i32>::new();
        eqx.orbit_reps(&mut reps);

        // Join onesies where possible.  This should probably be more efficient.

        for (&e1, info1) in exacts.iter().zip(to_infos.iter()).take(nexacts) {
            let ex1 = &exact_clonotypes[e1];
            if ex1.share.solo() {
                let mut is = Vec::<usize>::new();
                for (&e2, info2) in exacts.iter().take(nexacts).zip(to_infos.iter()) {
                    let ex2 = &exact_clonotypes[e2];
                    if ex2.share.solo() {
                        if ex1.share[0].seq == ex2.share[0].seq {
                            eqx.join(info1[0] as i32, info2[0] as i32);
                        }
                    } else {
                        for j in 0..ex2.share.len() {
                            if ex2.share[j].seq == ex1.share[0].seq {
                                is.push(info2[0]);
                            }
                        }
                    }
                }
                let mut rs = Vec::<usize>::new();
                for &ij in &is {
                    rs.push(eqx.class_id(ij as i32) as usize);
                }
                unique_sort(&mut rs);
                if rs.solo() {
                    eqx.join(info1[0] as i32, is[0] as i32);
                }
            }
        }

        // Divide the orbit if needed.

        if eqx.norbits() == 1 {
            orbits2.push(o.clone());
        } else {
            let mut repsx = Vec::<i32>::new();
            eqx.orbit_reps(&mut repsx);
            for r in repsx {
                let mut ox = Vec::<i32>::new();
                eqx.orbit(r, &mut ox);
                let mut o2 = Vec::<i32>::new();
                for ko in ox {
                    o2.push(o[ko as usize]);
                }
                orbits2.push(o2);
            }
        }
    }
    *orbits = orbits2;
}
