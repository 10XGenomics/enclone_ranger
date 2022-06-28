// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::print_utils1::color_codon;
use bio_edit::alignment::pairwise::Aligner;
use bio_edit::alignment::AlignmentOperation::{Del, Ins, Match, Subst};
use enclone_core::allowed_vars::{CVARS_ALLOWED, CVARS_ALLOWED_PCELL};
use enclone_core::defs::{ColInfo, EncloneControl, ExactClonotype};
use enclone_proto::types::DonorReferenceItem;
use io_utils::fwriteln;
use itertools::Itertools;
use std::fmt::Write as _;
use std::io::Write;
use string_utils::strme;
use vdj_ann::refx::RefData;
use vector_utils::{bin_member, erase_if, make_freq, next_diff, sort_sync2, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn comp_edit(
    ex: &ExactClonotype,
    mid: usize,
    col: usize,
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    rsi: &ColInfo,
) -> (usize, String) {
    let mut comp = 1000000;
    let mut edit = String::new();
    let td = &ex.share[mid];
    let tig = &td.seq;
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
    let mut aligner = Aligner::new(-6, -1, &score);

    // Go through passes.  If IGH/TRB, we go through every D segment.  Otherwise
    // there is just one pass.

    let mut z = 1;
    if ex.share[mid].left {
        z = refdata.ds.len();
    }
    let mut ds = Vec::<usize>::new();
    let mut counts = Vec::<usize>::new();
    for di in 0..z {
        let mut d = 0;
        if ex.share[mid].left {
            d = refdata.ds[di];
        }

        // Start to build reference concatenation.  First append the V segment.

        let mut concat = Vec::<u8>::new();
        let mut vref = refdata.refs[rsi.vids[col]].to_ascii_vec();
        if rsi.vpids[col].is_none() {
        } else {
            vref = dref[rsi.vpids[col].unwrap()].nt_sequence.clone();
        }
        concat.append(&mut vref.clone());

        // Append the D segment if IGH/TRB.

        if ex.share[mid].left {
            let mut x = refdata.refs[d].to_ascii_vec();
            concat.append(&mut x);
        }

        // Append the J segment.

        let mut x = refdata.refs[rsi.jids[col]].to_ascii_vec();
        concat.append(&mut x);

        // Align the V..J sequence on the contig to the reference concatenation.

        let al = aligner.semiglobal(tig, &concat);
        let mut m = 0;
        let mut pos = al.xstart;
        let mut rpos = (al.ystart as isize) - (vref.len() as isize);
        let mut count = 0;
        let start = td.cdr3_start - td.ins_len();
        let stop = td.j_stop - td.v_start;
        let mut edits = Vec::<String>::new();
        while m < al.operations.len() {
            let n = next_diff(&al.operations, m);
            match al.operations[m] {
                Match => {
                    pos += 1;
                    rpos += 1;
                }
                Subst => {
                    if pos >= start && pos < stop {
                        count += 1;
                        edits.push(format!("S{}", rpos));
                    }
                    pos += 1;
                    rpos += 1;
                }
                Del => {
                    if pos >= start && pos < stop {
                        count += 1;
                        edits.push(format!("D{}:{}", rpos, n - m));
                    }
                    pos += n - m;
                    m = n - 1;
                }
                Ins => {
                    if pos >= start && pos < stop {
                        count += 1;
                        edits.push(format!("I{}:{}", rpos, n - m));
                    }
                    rpos += (n - m) as isize;
                    m = n - 1;
                }
                _ => {}
            };
            m += 1;
        }
        counts.push(count);
        ds.push(d);
        if count < comp {
            comp = count;
            edit = format!("{}", edits.iter().format("•"));
        }
    }
    sort_sync2(&mut counts, &mut ds);
    let mut comp = 0;
    if !counts.is_empty() {
        comp = counts[0];
    }
    (comp, edit)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Form the consensus CDR3 DNA sequence for a clonotype.  We define this by finding, for each
// chain, and each codon in CDR3, the most frequent one.  These codons are then chained together
// to form a DNA sequence.  If no codon has > 50% of the total, we report a tie, as XXX.
// This code is only used in processing COLOR=codon-diffs.

pub fn consensus_codon_cdr3(
    rsi: &ColInfo,
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
) -> Vec<Vec<u8>> {
    let mut cons = Vec::<Vec<u8>>::new();
    let nexacts = exacts.len();
    let cols = rsi.vids.len();
    for cx in 0..cols {
        let mut cdr3s = Vec::<Vec<u8>>::new();
        for u in 0..nexacts {
            if let Some(m) = rsi.mat[cx][u] {
                cdr3s.push(
                    exact_clonotypes[exacts[u]].share[m]
                        .cdr3_dna
                        .as_bytes()
                        .to_vec(),
                );
            }
        }
        let n = cdr3s[0].len();
        let mut con = Vec::<u8>::new();
        for i in (0..n).step_by(3) {
            let mut codons = Vec::<Vec<u8>>::new();
            for cdr3 in &cdr3s {
                codons.push(cdr3[i..i + 3].to_vec());
            }
            codons.sort();
            let mut freq = Vec::<(u32, Vec<u8>)>::new();
            make_freq(&codons, &mut freq);
            if freq[0].0 as usize * 2 > codons.len() {
                con.append(&mut freq[0].1.clone());
            } else {
                con.append(&mut b"XXX".to_vec());
            }
        }
        cons.push(con);
    }
    cons
}

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
            let mut chain_descrip = format!("{}|{}", vdescrip, refdata.name[vid]);
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

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Add header text to mlog.

pub fn add_header_text(
    ctl: &EncloneControl,
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    rord: &[usize],
    mat: &[Vec<Option<usize>>],
    mut mlog: &mut Vec<u8>,
) {
    let nexacts = exacts.len();
    let cols = mat.len();
    for (cx, mcx) in mat.iter().take(cols).enumerate() {
        let (mut vref, mut jref) = (Vec::<u8>::new(), Vec::<u8>::new());
        for (&m, &e) in mcx.iter().zip(exacts.iter()).take(nexacts) {
            if let Some(m) = m {
                vref = exact_clonotypes[e].share[m].vs.to_ascii_vec();
                jref = exact_clonotypes[e].share[m].js.to_ascii_vec();
            }
        }
        let mut seqs = Vec::<Vec<u8>>::new();
        let mut full_seqs = Vec::<Vec<u8>>::new();
        for u in 0..nexacts {
            let ex = &exact_clonotypes[exacts[rord[u]]];
            if let Some(m) = mat[cx][rord[u]] {
                seqs.push(ex.share[m].seq_del.clone());
                full_seqs.push(ex.share[m].full_seq.clone());
            } else {
                seqs.push(Vec::<u8>::new());
                full_seqs.push(Vec::<u8>::new());
            }
        }
        let mut simple = false;
        if ctl.clono_print_opt.note_simple && vref.len() + jref.len() >= seqs[0].len() {
            let n = seqs[0].len() - jref.len();
            let mut vj = vref[0..n].to_vec();
            vj.append(&mut jref.clone());
            if vj == seqs[0] {
                simple = true;
            }
        }
        if ctl.clono_print_opt.seqc || ctl.clono_print_opt.full_seqc || simple {
            fwriteln!(&mut mlog, "CHAIN {}", cx + 1);
        }
        if ctl.clono_print_opt.seqc {
            fwriteln!(&mut mlog, "• {}", strme(&seqs[0]));
        }
        if ctl.clono_print_opt.full_seqc {
            fwriteln!(&mut mlog, "• {}", strme(&full_seqs[0]));
        }
        if simple {
            fwriteln!(&mut mlog, "• This chain is simple.");
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Insert universal and donor reference rows.
// Possibly buggy WRT reference indels that we insert.

pub fn insert_reference_rows(
    ctl: &EncloneControl,
    rsi: &ColInfo,
    show_aa: &[Vec<usize>],
    field_types: &[Vec<u8>],
    refdata: &RefData,
    dref: &[DonorReferenceItem],
    row1: &[String],
    drows: &mut Vec<Vec<String>>,
    rows: &mut Vec<Vec<String>>,
    exacts: &[usize],
    exact_clonotypes: &[ExactClonotype],
    peer_groups: &[Vec<(usize, u8, u32)>],
    cdr3_con: &[Vec<u8>],
) {
    let cols = rsi.seq_del_lens.len();
    if !drows.is_empty() {
        for pass in 1..=2 {
            let mut row = Vec::<String>::new();
            if pass == 1 {
                row.push("reference".to_string());
            } else {
                row.push("donor ref".to_string());
            }
            for _ in 1..row1.len() {
                row.push("\\ext".to_string());
            }
            for cz in 0..cols {
                let mut refseq = Vec::<u8>::new();
                let mut vlen: usize;
                let vseq: Vec<u8>;
                if pass == 1 || rsi.vpids[cz].is_none() {
                    vlen = refdata.refs[rsi.vids[cz]].len();
                    vseq = refdata.refs[rsi.vids[cz]].to_ascii_vec();
                } else {
                    vlen = dref[rsi.vpids[cz].unwrap()].nt_sequence.len();
                    vseq = dref[rsi.vpids[cz].unwrap()].nt_sequence.clone();
                }
                let mut trim = ctl.heur.ref_v_trim;
                vlen -= trim;
                let mut jlen = refdata.refs[rsi.jids[cz]].len() - trim;
                let jseq = refdata.refs[rsi.jids[cz]].to_ascii_vec();
                let mut gap = rsi.seq_del_lens[cz] as isize - vlen as isize - jlen as isize;

                if gap < -2 * (trim as isize) {
                    let mut bcs = Vec::<String>::new();
                    for u in 0..exacts.len() {
                        let ex = &exact_clonotypes[exacts[u]];
                        for i in 0..ex.clones.len() {
                            bcs.push(ex.clones[i][0].barcode.clone());
                        }
                    }
                    bcs.sort();
                    eprintln!("\ncz = {}", cz);
                    eprintln!("pass = {}", pass);
                    eprintln!("seq_del.len() = {}", rsi.seq_del_lens[cz]);
                    eprintln!("vlen = {}", vlen);
                    eprintln!("jlen = {}", jlen);
                    eprintln!("gap = seq_del.len() - vlen - jlen");
                    panic!(
                        "Something is wrong because gap is {}, which is negative.\n\
                        This is happening for the clonotype with these barcodes:\n{}.",
                        gap,
                        bcs.iter().format(",")
                    );
                }

                if gap < 0 {
                    let mut ptrim = (-gap) / 2;
                    if (-gap) % 2 == 1 {
                        ptrim += 1;
                    }
                    vlen += ptrim as usize;
                    jlen += ptrim as usize;
                    gap += 2 * ptrim;
                    trim -= ptrim as usize;
                }

                refseq.extend(&vseq[..vlen]);
                let gap = gap as usize;
                refseq.resize(refseq.len() + gap, b'-');
                refseq.extend(&jseq[trim..trim + jlen]);
                let mut refx = String::new();
                let mut last_color = "black".to_string();
                for k in 0..show_aa[cz].len() {
                    let p = show_aa[cz][k];
                    if k > 0
                        && field_types[cz][k] != field_types[cz][k - 1]
                        && !ctl.gen_opt.nospaces
                    {
                        refx += " ";
                    }
                    if 3 * p + 3 > refseq.len() || refseq[3 * p..3 * p + 3].contains(&b'-') {
                        refx += "◦";
                    } else {
                        let x = &peer_groups[rsi.vids[cz]];
                        let last = k == show_aa[cz].len() - 1;
                        let log = color_codon(
                            ctl,
                            &refseq,
                            &Vec::new(),
                            x,
                            cz,
                            0,
                            p,
                            0,
                            &mut last_color,
                            last,
                            cdr3_con,
                            exacts,
                            exact_clonotypes,
                        );
                        refx += strme(&log);
                    }
                }
                row.push(refx);
                for _ in 1..rsi.cvars[cz].len() {
                    row.push("".to_string());
                }
            }
            rows.push(row);
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Process COMPLETE.

pub fn process_complete(
    ctl: &EncloneControl,
    nexacts: usize,
    bads: &mut [bool],
    mat: &[Vec<Option<usize>>],
) {
    let cols = mat.len();
    if ctl.gen_opt.complete {
        let mut used = vec![false; cols];
        for (u, &b) in bads.iter().take(nexacts).enumerate() {
            if !b {
                for (used, m) in used.iter_mut().zip(mat.iter()).take(cols) {
                    if m[u].is_some() {
                        *used = true;
                    }
                }
            }
        }
        for (&used, mat) in used.iter().zip(mat.iter()).take(cols) {
            if used {
                for (b, m) in bads.iter_mut().take(nexacts).zip(mat.iter()) {
                    if m.is_none() {
                        *b = true;
                    }
                }
            }
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Identify certain extra parseable variables.  These arise from parameterizable cvars.

pub fn get_extra_parseables<'a>(ctl: &'a EncloneControl, pcols_sort: &'a [String]) -> Vec<&'a str> {
    let mut extra_parseables = Vec::<&str>::new();
    let mut exclusions = ctl
        .clono_print_opt
        .cvars
        .iter()
        .map(String::as_str)
        .collect::<Vec<_>>();
    for v in CVARS_ALLOWED {
        exclusions.push(v);
    }
    for v in CVARS_ALLOWED_PCELL {
        exclusions.push(v);
    }
    unique_sort(&mut exclusions);
    for x in pcols_sort.iter() {
        let chars = x.char_indices().collect::<Vec<_>>();
        let mut trim = 0;
        for c in chars.iter().rev() {
            if !c.1.is_digit(10) {
                break;
            }
            trim += 1;
        }
        if trim > 0 {
            let v = &x[..chars[chars.len() - trim - 1].0];
            if !bin_member(&exclusions, &v) {
                extra_parseables.push(v);
            }
        }
    }
    unique_sort(&mut extra_parseables);
    extra_parseables
}
