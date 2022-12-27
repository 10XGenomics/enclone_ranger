// Copyright (c) 2022 10X Genomics, Inc. All rights reserved.

// Analyze donor reference and exit.
//
// This displays tables, which are somewhat mangled unless there are at least four donors.

use debruijn::dna_string::DnaString;
use enclone_core::defs::EncloneControl;
use itertools::Itertools;
use std::cmp::min;
use std::env;
use string_utils::TextUtils;
use tables::print_tabular_vbox;
use vdj_ann::refx::{make_vdj_ref_data_core, RefData};
use vector_utils::{bin_position, erase_if, next_diff1_2, next_diff1_3, unique_sort};

pub fn analyze_donor_ref(
    refdata: &RefData,
    ctl: &EncloneControl,
    alt_refs: &[(usize, usize, DnaString, usize, bool)],
) {
    // Analyze donor reference.

    if !ctl.gen_opt.external_ref.is_empty() {
        if ctl.gen_opt.echo {
            let args: Vec<String> = env::args().collect();
            println!("{}", args.iter().format(" "));
        }
        let mut erefdata = RefData::new();
        let f = std::fs::read_to_string(&ctl.gen_opt.external_ref).unwrap();
        make_vdj_ref_data_core(&mut erefdata, &f, "", true, true, None);
        let mut refs = Vec::<(String, String, Vec<u8>)>::new(); // {(gene, allele, seq)}

        // Store the external (IMGT) alleles.

        for i in 0..erefdata.refs.len() {
            if erefdata.is_v(i) {
                let allele = erefdata.rheaders_orig[i].between("*", " ");
                refs.push((
                    erefdata.name[i].clone(),
                    allele.to_string(),
                    erefdata.refs[i].to_ascii_vec(),
                ));
            }
        }

        // Store the enclone reference alleles.

        for i in 0..refdata.refs.len() {
            if refdata.is_v(i) {
                refs.push((
                    refdata.name[i].clone(),
                    format!("uref{i}"),
                    refdata.refs[i].to_ascii_vec(),
                ));
            }
        }

        // Store the donor reference alleles;

        for (i, ar) in alt_refs.iter().enumerate() {
            let donor = ar.0;
            let ref_id = ar.1;
            let name = &refdata.name[ref_id];
            let alt_seq = &ar.2;
            refs.push((
                name.clone(),
                format!("dref{i}_{donor}"),
                alt_seq.to_ascii_vec(),
            ));
        }

        // Sort the alleles and group by gene.

        refs.sort();
        let mut i = 0;
        while i < refs.len() {
            let j = next_diff1_3(&refs, i as i32) as usize;
            let gene = &refs[i].0;
            let mut alleles = Vec::<(&[u8], &str)>::new(); // (sequence, name)
            let mut have_alt = false;
            for r in &refs[i..j] {
                if r.1.starts_with("dref") {
                    have_alt = true;
                }
                alleles.push((r.2.as_ref(), r.1.as_str()));
            }

            // Delete reference alleles having very low count relative to others.

            let mut to_delete = vec![false; alleles.len()];
            let mut mm = 0;
            for ak in &alleles {
                if ak.1.starts_with("dref") {
                    let ii = ak.1.between("dref", "_").force_usize();
                    mm = std::cmp::max(mm, alt_refs[ii].3);
                }
            }
            for (ak, d) in alleles.iter().zip(to_delete.iter_mut()) {
                if ak.1.starts_with("dref") {
                    let ii = ak.1.between("dref", "_").force_usize();
                    if alt_refs[ii].4 && alt_refs[ii].3 * 10 < mm {
                        *d = true;
                    }
                }
            }
            erase_if(&mut alleles, &to_delete);

            // Proceed.

            if have_alt {
                // Truncate alleles so that they all have the same length.

                let mut m = 1000000;
                for ar in &alleles {
                    m = min(m, ar.0.len());
                }
                for ar in alleles.iter_mut() {
                    if ar.0.len() > m {
                        ar.0 = &ar.0[..m];
                    }
                }

                // Now alleles = all the alleles for one gene, and there is at least one
                // donor reference allele.  Combine identical alleles, and reorder.

                alleles.sort();
                let mut allelesg = Vec::<(Vec<&str>, &[u8])>::new();
                let mut r = 0;
                while r < alleles.len() {
                    let s = next_diff1_2(&alleles, r as i32) as usize;
                    let mut names = Vec::<&str>::new();
                    for at in &alleles[r..s] {
                        names.push(at.1);
                    }
                    allelesg.push((names, alleles[r].0));
                    r = s;
                }

                // Find the positions at which the alleles differ, and make variant matrix.

                let mut dp = Vec::<usize>::new();
                for p in 0..m {
                    let mut bases = Vec::<u8>::new();
                    for ar in &allelesg {
                        bases.push(ar.1[p]);
                    }
                    unique_sort(&mut bases);
                    if bases.len() > 1 {
                        dp.push(p);
                    }
                }
                let mut dm = vec![vec![0_u8; dp.len()]; allelesg.len()];
                for (ar, dmr) in allelesg.iter().zip(dm.iter_mut()) {
                    for (dmu, &dpu) in dmr.iter_mut().zip(dp.iter()) {
                        *dmu = ar.1[dpu];
                    }
                }

                // Make donor matrix.

                let ndonors = ctl.origin_info.donor_list.len();
                let mut dd = vec![vec![false; ndonors]; allelesg.len()];
                for r in 0..allelesg.len() {
                    for k in 0..allelesg[r].0.len() {
                        let n = &allelesg[r].0[k];
                        if n.starts_with("dref") {
                            let d = n.after("_").force_usize();
                            dd[r][d] = true;
                        }
                    }
                }

                // Make IMGT matrix.

                let mut imgts = Vec::<&str>::new();
                for ar in &allelesg {
                    for n in ar.0.iter() {
                        if !n.starts_with('d') && !n.starts_with('u') {
                            imgts.push(n);
                        }
                    }
                }
                unique_sort(&mut imgts);
                let nimgt = imgts.len();
                let mut im = vec![vec![false; nimgt]; allelesg.len()];
                for ar in &allelesg {
                    for n in &ar.0 {
                        let p = bin_position(&imgts, n);
                        if p >= 0 {
                            im[r][p as usize] = true;
                        }
                    }
                }

                // Make table, if it won't be too wide.

                let mut log = String::new();
                if dp.len() <= 20 {
                    let mut rows = vec![
                        {
                            let mut row = vec!["allele".to_string(), "donor".to_string()];
                            for _ in 0..ndonors - 1 {
                                row.push("\\ext".to_string());
                            }
                            if nimgt > 0 {
                                row.push("IMGT".to_string());
                                for _ in 0..nimgt - 1 {
                                    row.push("\\ext".to_string());
                                }
                            }
                            if !dp.is_empty() {
                                row.push("position".to_string());
                                for _ in 0..dp.len() - 1 {
                                    row.push("\\ext".to_string());
                                }
                            }
                            row
                        },
                        {
                            let mut row = vec!["".to_string()];
                            row.append(&mut vec![
                                "\\hline".to_string();
                                ndonors + nimgt + dp.len()
                            ]);
                            row
                        },
                        {
                            let mut row = vec!["".to_string()];
                            for d in 0..ndonors {
                                row.push(format!("{}", d + 1));
                            }
                            for im in imgts {
                                row.push(im.to_string());
                            }
                            for &u in &dp {
                                row.push(u.to_string());
                            }
                            row
                        },
                        vec!["\\hline".to_string(); ndonors + nimgt + dp.len() + 1],
                    ];
                    for (r, alleleg) in allelesg.into_iter().enumerate() {
                        let mut row = Vec::<String>::new();
                        let allele_name = (b'A' + r as u8) as char;
                        let mut an = String::new();
                        an.push(allele_name);
                        for n in alleleg.0.iter() {
                            if n.starts_with("uref") {
                                an.push('*');
                                break;
                            }
                        }
                        row.push(an);
                        for d in 0..ndonors {
                            if dd[r][d] {
                                row.push("▓".to_string());
                            } else {
                                row.push(" ".to_string());
                            }
                        }
                        for k in 0..nimgt {
                            if im[r][k] {
                                row.push("▓".to_string());
                            } else {
                                row.push(" ".to_string());
                            }
                        }
                        for &u in &dp {
                            row.push((alleleg.1[u] as char).to_string());
                        }
                        rows.push(row);
                    }
                    let mut just = b"l|".to_vec();
                    just.extend(vec![b'l'; ndonors]);
                    if nimgt > 0 {
                        just.push(b'|');
                        just.extend(vec![b'l'; nimgt]);
                    }
                    if !dp.is_empty() {
                        just.push(b'|');
                        just.extend(vec![b'l'; dp.len()]);
                    }
                    print_tabular_vbox(&mut log, &rows, 1, &just, false, false);
                }

                // Print.

                println!("\nworking on {gene}, have {} seqs", alleles.len());
                println!(
                    "alleles differ at {} positions = {}",
                    dp.len(),
                    dp.iter().format(",")
                );
                if !log.is_empty() {
                    log.truncate(log.len() - 1);
                    println!("\n{log}");
                    println!("* = a universal reference\n");
                }
                for m1 in 0..alleles.len() {
                    for m2 in m1 + 1..alleles.len() {
                        let a1 = &alleles[m1];
                        let a2 = &alleles[m2];
                        let mut diffs = 0;
                        for p in 0..min(a1.0.len(), a2.0.len()) {
                            if a1.0[p] != a2.0[p] {
                                diffs += 1;
                            }
                        }
                        println!(
                            "{} = {} vs {} = {} ==> {} diffs",
                            m1 + 1,
                            a1.1,
                            m2 + 1,
                            a2.1,
                            diffs
                        );
                    }
                }
                for a1 in &alleles {
                    let mut best = 1_000_000;
                    if !a1.1.starts_with("dref") {
                        continue;
                    }
                    for a2 in &alleles {
                        if a2.1.starts_with("dref") {
                            continue;
                        }
                        let mut diffs = 0;
                        for p in 0..min(a1.0.len(), a2.0.len()) {
                            if a1.0[p] != a2.0[p] {
                                diffs += 1;
                            }
                        }
                        best = min(best, diffs);
                    }
                    println!("{} is distance {best} from a reference", a1.1);
                }
            }
            i = j;
        }
        std::process::exit(0);
    }
}
