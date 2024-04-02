// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use std::fmt::Write as _;
use vdj_ann::refx::RefData;
use vector_utils::{erase_if, make_freq};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Define reference segment identifiers, one per chain.  For the V segment, there is
// also an optional donor reference sequence alignment.  (For now.  We might do the
// same thing later for J segments.)
// ◼ 1. We vote using the number of cells, whereas a better way would be to test all
// ◼    the alternatives to find the best match.
// ◼ 2. Defining a constant region identifier for the entire clonotype is
// ◼    biologically dubious.
// ◼ 3. Maybe we only need to do this for pass 2.

pub fn define_column_info(
    ctl: &EncloneControl,
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    mat: &[Vec<Option<usize>>],
    refdata: &RefData,
) -> ColInfo {
    let cols = mat.len();

    // Define cvars.

    let mut cvars = Vec::<Vec<String>>::new();
    for m in mat.iter().take(cols) {
        let mut have_notes = false;
        let mut left = false;
        for (&e, &m) in exacts.iter().zip(m.iter()) {
            let ex = &exact_clonotypes[e];
            if let Some(m) = m {
                let ex = &ex.share[m];
                if ex.left {
                    left = true;
                }
                if !ex.vs_notesx.is_empty() {
                    have_notes = true;
                }
            }
        }
        let mut cv = Vec::<String>::new();
        for i in 0..ctl.clono_print_opt.cvars.len() {
            let var = &ctl.clono_print_opt.cvars[i];
            if var == "notes" && !have_notes {
                continue;
            }
            if !left
                && (var == "d1_name"
                    || var == "d2_name"
                    || var == "d_delta"
                    || var == "d_Δ"
                    || var == "d1_score"
                    || var == "d2_score")
            {
                continue;
            }
            cv.push(var.to_string());
        }
        cvars.push(cv);
    }

    // Compute CDR3 starts, etc.

    let mut fr1_starts = Vec::<usize>::new();
    let mut fr2_starts = Vec::<Option<usize>>::new();
    let mut fr3_starts = Vec::<Option<usize>>::new();
    let mut cdr1_starts = Vec::<Option<usize>>::new();
    let mut cdr2_starts = Vec::<Option<usize>>::new();
    let mut cdr3_starts = Vec::<usize>::new();
    let mut cdr3_lens = Vec::<usize>::new();
    let mut seq_lens = Vec::<usize>::new();
    let mut seq_del_lens = Vec::<usize>::new();
    for m in mat.iter().take(cols) {
        for (&e, &m) in exacts.iter().zip(m.iter()) {
            let ex = &exact_clonotypes[e];
            if let Some(m) = m {
                let exm = &ex.share[m];
                cdr3_lens.push(exm.cdr3_aa.len());
                seq_lens.push(exm.seq.len());
                seq_del_lens.push(exm.seq_del.len());

                // The logic below with testing i < start while incrementing start seems fishy.

                if let Some(mut start) = exm.cdr1_start {
                    for (i, c) in exm.seq_del.iter().enumerate() {
                        if i < start && *c == b'-' {
                            start += 1;
                        }
                    }
                    cdr1_starts.push(Some(start));
                } else {
                    cdr1_starts.push(None);
                }
                if let Some(mut start) = exm.cdr2_start {
                    for (i, c) in exm.seq_del.iter().enumerate() {
                        if i < start && *c == b'-' {
                            start += 1;
                        }
                    }
                    cdr2_starts.push(Some(start));
                } else {
                    cdr2_starts.push(None);
                }
                let mut start = exm.fr1_start;
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                fr1_starts.push(start);
                if exm.fr2_start.is_some() {
                    let mut start = exm.fr2_start.unwrap();
                    for (i, c) in exm.seq_del.iter().enumerate() {
                        if i < start && *c == b'-' {
                            start += 1;
                        }
                    }
                    fr2_starts.push(Some(start));
                } else {
                    fr2_starts.push(None);
                }
                if exm.fr3_start.is_some() {
                    let mut start = exm.fr3_start.unwrap();
                    for (i, c) in exm.seq_del.iter().enumerate() {
                        if i < start && *c == b'-' {
                            start += 1;
                        }
                    }
                    fr3_starts.push(Some(start));
                } else {
                    fr3_starts.push(None);
                }
                let mut start = exm.cdr3_start - exm.ins_len();
                for (i, c) in exm.seq_del.iter().enumerate() {
                    if i < start && *c == b'-' {
                        start += 1;
                    }
                }
                cdr3_starts.push(start);
                break;
            }
        }
    }

    // Compute reference info.

    let mut uids = vec![None; cols];
    let mut vids = vec![0; cols];
    let mut vpids = vec![None; cols];
    let mut vpids_d = vec![None; cols];
    let mut vpids_a = vec![None; cols];
    let mut dids = vec![None; cols];
    let mut jids = vec![0; cols];
    let mut cids = vec![None; cols];
    let mut left = vec![false; cols];
    for (m, (left, ((uids, (vids, (vpids, (vpids_d, vpids_a)))), (dids, (jids, cids))))) in mat
        .iter()
        .zip(
            left.iter_mut().zip(
                uids.iter_mut()
                    .zip(
                        vids.iter_mut().zip(
                            vpids
                                .iter_mut()
                                .zip(vpids_d.iter_mut().zip(vpids_a.iter_mut())),
                        ),
                    )
                    .zip(dids.iter_mut().zip(jids.iter_mut().zip(cids.iter_mut()))),
            ),
        )
        .take(cols)
    {
        let mut u = Vec::<usize>::new();
        let mut v = Vec::<usize>::new();
        let mut vp = Vec::<(usize, Option<usize>, Option<usize>, Option<usize>)>::new();
        let mut d = Vec::<usize>::new();
        let mut j = Vec::<usize>::new();
        let mut c = Vec::<usize>::new();
        for (&clonotype_id, &m) in exacts.iter().zip(m.iter()) {
            let ex = &exact_clonotypes[clonotype_id];
            if let Some(m) = m {
                let x = &ex.share[m];
                if x.left {
                    *left = true;
                }
                let ncells = ex.ncells();
                if let Some(u_ref_id) = x.u_ref_id {
                    u.resize(u.len() + ncells, u_ref_id);
                }
                // This is not actually correct.  It copies the consensus V gene assignment
                // for an exact subclonotype, rather than fetch the per cell entries.  However
                // it would be very rare for this to make a difference.
                v.resize(v.len() + ncells, x.v_ref_id);
                vp.resize(
                    vp.len() + ncells,
                    (
                        x.v_ref_id,
                        x.v_ref_id_donor,
                        x.v_ref_id_donor_donor,
                        x.v_ref_id_donor_alt_id,
                    ),
                );
                if let Some(d_ref_id) = x.d_ref_id {
                    d.resize(d.len() + ncells, d_ref_id);
                }
                j.resize(j.len() + ncells, x.j_ref_id);
                if let Some(c_ref_id) = x.c_ref_id {
                    c.resize(c.len() + ncells, c_ref_id);
                }
            }
        }
        u.sort_unstable();
        v.sort_unstable();
        vp.sort();
        d.sort_unstable();
        j.sort_unstable();
        c.sort_unstable();
        let mut uf = Vec::<(u32, usize)>::new();
        make_freq(&u, &mut uf);
        if !uf.is_empty() {
            *uids = Some(uf[0].1);
        }
        let mut vf = Vec::<(u32, usize)>::new();
        make_freq(&v, &mut vf);
        *vids = vf[0].1;
        let mut to_delete = vec![false; vp.len()];
        for i in 0..vp.len() {
            if vp[i].0 != *vids {
                to_delete[i] = true;
            }
        }
        erase_if(&mut vp, &to_delete);
        let mut vpf = Vec::<(u32, (usize, Option<usize>, Option<usize>, Option<usize>))>::new();
        make_freq(&vp, &mut vpf);
        *vpids = (vpf[0].1).1;
        *vpids_d = (vpf[0].1).2;
        *vpids_a = (vpf[0].1).3;
        let mut df = Vec::<(u32, usize)>::new();
        make_freq(&d, &mut df);
        if !df.is_empty() {
            *dids = Some(df[0].1);
        }
        let mut jf = Vec::<(u32, usize)>::new();
        make_freq(&j, &mut jf);
        *jids = jf[0].1;
        let mut cf = Vec::<(u32, usize)>::new();
        make_freq(&c, &mut cf);
        if !cf.is_empty() {
            *cids = Some(cf[0].1);
        }
    }

    // Compute seqss and seqss_amino.

    let mut seqss = Vec::<Vec<Vec<u8>>>::new();
    let mut seqss_amino = Vec::<Vec<Vec<u8>>>::new();
    let nexacts = exacts.len();
    for cx in mat.iter().take(cols) {
        let mut seqs = Vec::<Vec<u8>>::new();
        let mut seqs_amino = Vec::<Vec<u8>>::new();
        for (&m, &e) in cx.iter().zip(exacts.iter()).take(nexacts) {
            if let Some(m) = m {
                seqs.push(exact_clonotypes[e].share[m].seq_del.clone());
                seqs_amino.push(exact_clonotypes[e].share[m].seq_del_amino.clone());
            } else {
                seqs.push(Vec::<u8>::new());
                seqs_amino.push(Vec::<u8>::new());
            }
        }
        seqss.push(seqs.clone());
        seqss_amino.push(seqs_amino.clone());
    }

    // Show segment names.  We used ◼ as a separator character, but that does not render well
    // as a fixed-width character in Google Docs.  So we changed it to ◆.

    let chain_descrip = (0..cols)
        .map(|cx| {
            let vid = vids[cx];
            let mut vdescrip = format!("{}", refdata.id[vid]);
            if vpids[cx].is_some() {
                vdescrip = format!(
                    "{}.{}.{}",
                    vdescrip,
                    vpids_d[cx].unwrap() + 1,
                    vpids_a[cx].unwrap() + 1
                );
            }
            let mut chain_descrip = format!("{vdescrip}|{}", refdata.name[vid]);
            if let Some(did) = dids[cx] {
                write!(
                    chain_descrip,
                    " ◆ {}|{}",
                    refdata.id[did], refdata.name[did]
                )
                .unwrap();
            }
            let jid = jids[cx];
            write!(
                chain_descrip,
                " ◆ {}|{}",
                refdata.id[jid], refdata.name[jid]
            )
            .unwrap();
            chain_descrip
        })
        .collect();

    // Return.

    ColInfo {
        left,
        uids,
        vids,
        vpids,
        dids,
        jids,
        cids,
        fr1_starts,
        fr2_starts,
        fr3_starts,
        cdr1_starts,
        cdr2_starts,
        cdr3_starts,
        cdr3_lens,
        seq_lens,
        seq_del_lens,
        seqss,
        seqss_amino,
        chain_descrip,
        mat: Vec::<Vec<Option<usize>>>::new(),
        cvars,
    }
}
