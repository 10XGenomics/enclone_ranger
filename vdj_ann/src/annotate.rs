// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains code to annotate a contig, in the sense of finding alignments
// to VDJ reference contigs.  Also to find CDR3 sequences.  And some related things.

use crate::transcript::is_productive_contig;
use crate::{refx::RefData, transcript::ContigStatus};
use align_tools::affine_align;
use amino::{have_start, nucleotide_to_aminoacid_sequence};
use bio_edit::alignment::AlignmentOperation::{Del, Ins, Match, Subst, Xclip, Yclip};
use debruijn::{
    dna_string::{DnaString, DnaStringSlice},
    kmer::{Kmer12, Kmer20},
    Mer, Vmer,
};
use io_utils::{fwrite, fwriteln};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use stats_utils::percent_ratio;
use std::fmt::Write as _;
use std::{
    cmp::{max, min},
    fs::File,
    io::{BufWriter, Write},
};
use string_utils::{stringme, strme, TextUtils};
use vdj_types::{VdjChain, VdjRegion};
use vector_utils::{
    bin_member, erase_if, lower_bound1_3, next_diff12_4, next_diff1_2, next_diff1_3, next_diff1_5,
    reverse_sort, unique_sort, VecUtils,
};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// START CODONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn print_start_codon_positions(tig: &DnaString, log: &mut Vec<u8>) {
    let mut starts = Vec::<usize>::new();
    if tig.len() < 3 {
        return;
    }
    for i in 0..tig.len() - 3 {
        if have_start(tig, i) {
            starts.push(i);
        }
    }
    fwriteln!(log, "start codons at {}", starts.iter().format(", "));
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// ASSIGN CHAIN TYPE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Assign a chain type to a given DNA sequence b.
//
// The chain type is either -1, meaning unassigned, or an index into the vector
// "IGH","IGK","IGL","TRA","TRB","TRD","TRG"
// representing a forward alignment to the chain type, or 7 + such an index,
// representing a reverse alignment.
//
// This takes as input the sequence b, plus the following auxiliary data
// structures:
// * a 20-mer kmer lookup table for the VDJ reference sequences,
//   both TCR and BCR;
// * a classification vector that assigns each reference sequence to either a
//   chain type index or -1.
//
// ◼ Ns are incorrectly handled.  See lena 100349 for lots of examples.

pub fn chain_type(b: &DnaString, rkmers_plus_full_20: &[(Kmer20, i32, i32)], rtype: &[i32]) -> i8 {
    const N: usize = 7;
    let k = 20;
    if b.len() < k {
        return -1_i8;
    }
    let mut count_this = [0; 2 * N];
    let brc = b.rc();
    for l in 0..b.len() - k + 1 {
        let mut is_type = [false; 2 * N];
        for pass in 0..2 {
            let z = if pass == 1 { N } else { 0 };
            let x = if pass == 0 {
                b.get_kmer(l)
            } else {
                brc.get_kmer(l)
            };
            let low = lower_bound1_3(rkmers_plus_full_20, &x) as usize;
            for j in low..rkmers_plus_full_20.len() {
                if rkmers_plus_full_20[j].0 != x {
                    break;
                }
                let t = rkmers_plus_full_20[j].1 as usize;
                if rtype[t] >= 0 {
                    is_type[z + rtype[t] as usize] = true;
                }
            }
        }
        let mut nt = 0;
        for l in 0..2 * N {
            if is_type[l] {
                nt += 1;
            }
        }
        if nt == 1 {
            for l in 0..2 * N {
                if is_type[l] {
                    count_this[l] += 1;
                }
            }
        }
    }
    let m = *count_this.iter().max().unwrap();
    let best = count_this
        .iter()
        .enumerate()
        .rfind(|&v| *v.1 == m)
        .unwrap()
        .0;
    reverse_sort(&mut count_this);
    if count_this[0] > count_this[1] {
        best as i8
    } else {
        -1_i8
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// ANNOTATE SEQUENCES
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given a DnaString, enumerate matches to reference sequences.  Matches are
// defined to be gap-free alignments seeded on 12-mer matches, with mismatches
// allowed in the following cases:
//
// 1. Given two successive maximal perfect matches of length >= 12 having the same
//    offset, mismatches are allowed between them, so long as the error rate for
//    the extended match is at most 10%.
//
// 2. We always allow extension over a single mismatch so long as 5 perfectly
//    matching bases follow.
//
// However, we require a 20-mer match except for J regions.
// (see below for details)
//
// The structure of the output is:
// { ( start on sequence, match length, ref tig, start on ref tig, mismatches on sequence ) }.

pub fn annotate_seq(
    b: &DnaString,
    refdata: &RefData,
    ann: &mut Vec<(i32, i32, i32, i32, i32)>,
    allow_weak: bool,
    allow_improper: bool,
    abut: bool,
) {
    let mut log = Vec::<u8>::new();
    annotate_seq_core(
        b,
        refdata,
        ann,
        allow_weak,
        allow_improper,
        abut,
        &mut log,
        false,
    );
}

fn print_alignx(log: &mut Vec<u8>, a: &(i32, i32, i32, i32, Vec<i32>), refdata: &RefData) {
    let t = a.2 as usize;
    let l = a.0;
    let len = a.1;
    let p = a.3;
    let mis = a.4.len();
    fwriteln!(
        log,
        "{}-{} ==> {}-{} on {} (mis={})",
        l,
        l + len,
        p,
        p + len,
        refdata.rheaders[t],
        mis
    );
}

fn report_semis(
    verbose: bool,
    title: &str,
    semi: &[(i32, i32, i32, i32, Vec<i32>)],
    b_seq: &[u8],
    refs: &[DnaString],
    log: &mut Vec<u8>,
) {
    if verbose {
        fwriteln!(log, "\n{}\n", title);
        for s in semi {
            fwrite!(
                log,
                "t = {}, offset = {}, tig start = {}, ref start = {}, len = {}, mis = {}",
                s.0,
                s.1,
                s.2,
                s.1 + s.2,
                s.3,
                s.4.len(),
            );
            let t = s.0 as usize;
            let off = s.1;
            let tig_start = s.2;
            let ref_start = off + tig_start;
            let len = s.3;
            let mis = &s.4;
            let mut new_mis = Vec::<i32>::new();
            for j in 0..len {
                if b_seq[(tig_start + j) as usize] != refs[t].get((ref_start + j) as usize) {
                    new_mis.push(tig_start + j);
                }
            }
            if new_mis != *mis {
                fwriteln!(log, " [INVALID]");
                fwriteln!(log, "computed = {}", mis.iter().format(","));
                fwriteln!(log, "correct  = {}", new_mis.iter().format(","));
            } else {
                fwriteln!(log, "");
            }
        }
    }
}

pub fn annotate_seq_core(
    b: &DnaString,
    refdata: &RefData,
    ann: &mut Vec<(i32, i32, i32, i32, i32)>,
    allow_weak: bool,
    allow_improper: bool,
    abut: bool,
    log: &mut Vec<u8>,
    verbose: bool,
) {
    // The DNA string representation is inefficient because it stores bases as packed k-mers
    // which requires a lot of array bounds checks when unpacking which was a hot path
    // we found when profiling the CI job.  To avoid those in the inner
    // loop, we unpack it once here:
    let b_seq = b.to_bytes();
    // Unpack refdata.

    let refs = &refdata.refs;
    let rheaders = &refdata.rheaders;
    let rkmers_plus = &refdata.rkmers_plus;

    // Heuristic constants.

    const K: usize = 12;
    const MIN_PERF_EXT: usize = 5;
    const MAX_RATE: f64 = 0.15;

    // Find maximal perfect matches of length >= 20, or 12 for J regions, so long
    // as we have extension to a 20-mer with only one mismatch.
    //
    // perf = {(ref_id, ref_start - tig_start, tig_start, len)}

    let mut perf = Vec::<(i32, i32, i32, i32)>::new();
    if b.len() < K {
        return;
    }
    for l in 0..(b.len() - K + 1) {
        let x: Kmer12 = b.get_kmer(l);
        let low = lower_bound1_3(rkmers_plus, &x) as usize;
        for r in low..rkmers_plus.len() {
            if rkmers_plus[r].0 != x {
                break;
            }
            let t = rkmers_plus[r].1 as usize;
            let p = rkmers_plus[r].2 as usize;
            if l > 0 && p > 0 && b_seq[l - 1] == refs[t].get(p - 1) {
                continue;
            }
            let mut len = K;
            while l + len < b.len() && p + len < refs[t].len() {
                if b_seq[l + len] != refs[t].get(p + len) {
                    break;
                }
                len += 1;
            }
            let mut ok = len >= 20;
            if !ok && allow_weak {
                let mut ext1 = len + 1;
                let mut lx = l as i32 - 2;
                let mut px = p as i32 - 2;
                while lx >= 0 && px >= 0 {
                    if b_seq[lx as usize] != refs[t].get(px as usize) {
                        break;
                    }
                    ext1 += 1;
                    lx -= 1;
                    px -= 1;
                }
                let mut ext2 = len + 1;
                let mut lx = l + len + 1;
                let mut px = p + len + 1;
                while lx < b.len() && px < refs[t].len() {
                    if b_seq[lx] != refs[t].get(px) {
                        break;
                    }
                    ext2 += 1;
                    lx += 1;
                    px += 1;
                }
                if ext1 >= 20 || ext2 >= 20 {
                    ok = true;
                }
            }
            if ok {
                perf.push((t as i32, p as i32 - l as i32, l as i32, len as i32));
            }
        }
    }
    if verbose {
        fwriteln!(log, "\nINITIAL PERF ALIGNMENTS\n");
        for s in perf.iter() {
            fwriteln!(
                log,
                "t = {}, offset = {}, tig start = {}, ref start = {}, len = {}",
                s.0,
                s.1,
                s.2,
                s.1 + s.2,
                s.3,
            );
        }
    }

    // Find maximal perfect matches of length >= 10 that have the same offset as a perfect match
    // already found and are not equal to one of them.  But only do this if we already have at
    // least 150 bases aligned.

    let mut offsets = Vec::<(i32, i32)>::new();
    for i in 0..perf.len() {
        offsets.push((perf[i].0, perf[i].1));
    }
    unique_sort(&mut offsets);
    const MM_START: i32 = 150;
    const MM: i32 = 10;
    for (t, off) in offsets {
        let mut tig_starts = Vec::<i32>::new();
        let mut total = 0;
        for pi in &perf {
            if pi.0 == t && pi.1 == off {
                tig_starts.push(pi.2);
                total += pi.3;
            }
        }
        if total < MM_START {
            continue;
        }
        let r = refs[t as usize].to_bytes();
        let (mut l, mut p) = (0, off);
        while l <= b_seq.len() as i32 - MM {
            if p + MM > r.len() as i32 {
                break;
            }
            if p < 0 || b_seq[l as usize] != r[p as usize] {
                l += 1;
                p += 1;
            } else {
                let (mut lx, mut px) = (l + 1, p + 1);
                loop {
                    if lx >= b_seq.len() as i32 || px >= r.len() as i32 {
                        break;
                    }
                    if b_seq[lx as usize] != r[px as usize] {
                        break;
                    }
                    lx += 1;
                    px += 1;
                }
                let len = lx - l;
                if (MM..20).contains(&len) {
                    let mut known = false;
                    for k in 0..tig_starts.len() {
                        if l == tig_starts[k] {
                            known = true;
                            break;
                        }
                    }
                    if !known {
                        perf.push((t, p - l, l, len));
                    }
                }
                l = lx;
                p = px;
            }
        }
    }

    // Sort perfect matches.

    perf.sort_unstable();
    if verbose {
        fwriteln!(log, "\nPERF ALIGNMENTS\n");
        for s in perf.iter() {
            fwriteln!(
                log,
                "t = {}, offset = {}, tig start = {}, ref start = {}, len = {}",
                s.0,
                s.1,
                s.2,
                s.1 + s.2,
                s.3,
            );
        }
    }

    // Merge perfect matches.  We track the positions on b of mismatches.
    // semi = {(t, off, pos on b, len, positions on b of mismatches)}
    // where off = pos on ref - pos on b

    let mut semi = Vec::<(i32, i32, i32, i32, Vec<i32>)>::new();
    let mut i = 0;
    while i < perf.len() {
        let j = next_diff12_4(&perf, i as i32);
        let (t, off) = (perf[i].0, perf[i].1);
        let mut join = vec![false; j as usize - i];
        let mut mis = vec![Vec::<i32>::new(); j as usize - i];
        for k in i..j as usize - 1 {
            let (l1, len1) = (perf[k].2, perf[k].3);
            let (l2, len2) = (perf[k + 1].2, perf[k + 1].3);
            for z in l1 + len1..l2 {
                if b_seq[z as usize] != refs[t as usize].get((z + off) as usize) {
                    mis[k - i].push(z);
                }
            }

            // XXX:
            // println!( "\ntrying merge" );
            // printme!( t, l1, l2, len1, len2, mis[k-i].len() );

            if mis[k - i].len() as f64 / (l2 + len2 - l1) as f64 <= MAX_RATE {
                join[k - i] = true;
            }
        }
        let mut k1 = i;
        while k1 < j as usize {
            // let mut k2 = k1 + 1;
            let mut k2 = k1;
            let mut m = Vec::<i32>::new();
            // m.append( &mut mis[k1-i].clone() );
            while k2 < j as usize {
                // if !join[k2-i-1] { break; }
                if !join[k2 - i] {
                    break;
                }
                m.append(&mut mis[k2 - i].clone());
                k2 += 1;
            }
            semi.push((t, off, perf[k1].2, perf[k2].2 + perf[k2].3 - perf[k1].2, m));
            k1 = k2 + 1;
        }
        i = j as usize;
    }
    report_semis(verbose, "INITIAL SEMI ALIGNMENTS", &semi, &b_seq, refs, log);

    // Extend backwards and then forwards.

    for i in 0..semi.len() {
        let t = semi[i].0;
        let off = semi[i].1;
        let mut l = semi[i].2;
        let mut len = semi[i].3;
        let mut mis = semi[i].4.clone();
        while l > MIN_PERF_EXT as i32 && l + off > MIN_PERF_EXT as i32 {
            let mut ok = true;
            for j in 0..MIN_PERF_EXT {
                if b_seq[(l - j as i32 - 2) as usize]
                    != refs[t as usize].get((l + off - j as i32 - 2) as usize)
                {
                    ok = false;
                }
            }
            if !ok {
                break;
            }
            mis.push(l - 1);
            l -= MIN_PERF_EXT as i32 + 1;
            len += MIN_PERF_EXT as i32 + 1;
            while l > 0 && l + off > 0 {
                if b_seq[l as usize - 1] != refs[t as usize].get((l + off - 1) as usize) {
                    break;
                }
                l -= 1;
                len += 1;
            }
        }
        while l + len < (b.len() - MIN_PERF_EXT) as i32
            && l + len + off < (refs[t as usize].len() - MIN_PERF_EXT) as i32
        {
            let mut ok = true;
            for j in 0..MIN_PERF_EXT {
                if b_seq[(l + len + j as i32 + 1) as usize]
                    != refs[t as usize].get((l + off + len + j as i32 + 1) as usize)
                {
                    ok = false;
                }
            }
            if !ok {
                break;
            }
            mis.push(l + len);
            len += MIN_PERF_EXT as i32 + 1;
            while l + len < b.len() as i32 && l + off + len < refs[t as usize].len() as i32 {
                if b_seq[(l + len) as usize] != refs[t as usize].get((l + off + len) as usize) {
                    break;
                }
                len += 1;
            }
        }
        semi[i].2 = l;
        semi[i].3 = len;
        semi[i].4 = mis;
    }
    for i in 0..semi.len() {
        semi[i].4.sort_unstable();
    }

    // Add some 40-mers with the same offset having <= 6 mismatches.
    // semi = {(t, off, pos on b, len, positions on b of mismatches)}
    // where off = pos on ref - pos on b
    //
    // Note that implementation is asymmetric: we don't look to the right of p2, not for
    // any particularly good reason.
    //
    // This was added to get the heavy chain V segment of the mouse A20 cell line to be annotated.
    // This is dubious because the cell line is ~30 years old and of uncertain ancestry.  Thus
    // we're not sure if it arose from supported mouse strains or if the V segment might have
    // been corrupted during the growth of the cell line.  The A20 heavy chain V differs by 20%
    // from the reference.

    let mut i = 0;
    while i < semi.len() {
        let mut j = i + 1;
        let t = semi[i].0;
        let off = semi[i].1;
        while j < semi.len() {
            if semi[j].0 != t || semi[j].1 != off {
                break;
            }
            j += 1;
        }
        const L: i32 = 40;
        const MAX_DIFFS: usize = 6;
        let p1 = off + semi[i].2;
        // let p2 = off + semi[j-1].2 + semi[j-1].3;
        if -off >= 0 && p1 - off <= b_seq.len() as i32 {
            for p in 0..p1 - L {
                let l = p - off;
                let mut diffs = 0;
                for m in 0..L {
                    if b_seq[(l + m) as usize] != refs[t as usize].get((p + m) as usize) {
                        diffs += 1;
                        if diffs > MAX_DIFFS {
                            break;
                        }
                    }
                }
                if diffs <= MAX_DIFFS {
                    let mut x = Vec::<i32>::new();
                    for m in 0..L {
                        if b_seq[(l + m) as usize] != refs[t as usize].get((p + m) as usize) {
                            x.push(l + m);
                        }
                    }
                    semi.push((t, off, p - off, L, x));
                    break;
                }
            }
        }
        i = j;
    }
    semi.sort();

    // Allow extension over some mismatches on right if it gets us to the end on
    // the reference.  Ditto for left.
    // ◼ Not documented above.

    if allow_weak {
        let max_mis = 5;
        for i in 0..semi.len() {
            let t = semi[i].0;
            let off = semi[i].1;
            let l = semi[i].2;
            let mut len = semi[i].3;
            let mut mis = semi[i].4.clone();
            let mut mis_count = 0;
            while l + len < b_seq.len() as i32 && l + len + off < refs[t as usize].len() as i32 {
                if b_seq[(l + len) as usize] != refs[t as usize].get((l + off + len) as usize) {
                    mis.push(l + len);
                    mis_count += 1;
                }
                len += 1;
            }
            if mis_count <= max_mis && l + len + off == refs[t as usize].len() as i32 {
                semi[i].3 = len;
                semi[i].4 = mis;
            }
        }
        for i in 0..semi.len() {
            let t = semi[i].0;
            let off = semi[i].1;
            let mut l = semi[i].2;
            let mut len = semi[i].3;
            let mut mis = semi[i].4.clone();
            let mut mis_count = 0;
            while l > 0 && l + off > 0 {
                if b_seq[(l - 1_i32) as usize] != refs[t as usize].get((l + off - 1_i32) as usize) {
                    mis.push(l - 1);
                    mis_count += 1;
                }
                l -= 1;
                len += 1;
            }
            if mis_count <= max_mis && l + off == 0 {
                semi[i].2 = l;
                semi[i].3 = len;
                semi[i].4 = mis;
            }
        }
        for i in 0..semi.len() {
            semi[i].4.sort_unstable();
        }
    }
    report_semis(verbose, "SEMI ALIGNMENTS", &semi, &b_seq, refs, log);

    // Extend between match blocks.
    // ◼ This is pretty crappy.  What we should do instead is arrange the initial
    // ◼ extension between match blocks so it can be iterated.

    let mut to_delete = vec![false; semi.len()];
    for i1 in 0..semi.len() {
        let t1 = semi[i1].0;
        if t1 < 0 {
            continue;
        }
        let off1 = semi[i1].1;
        let (l1, len1) = (semi[i1].2, semi[i1].3);
        let mis1 = semi[i1].4.clone();
        for i2 in 0..semi.len() {
            let t2 = semi[i2].0;
            let off2 = semi[i2].1;
            if t2 != t1 || off2 != off1 {
                continue;
            }
            let (l2, len2) = (semi[i2].2, semi[i2].3);
            if l1 + len1 >= l2 {
                continue;
            }
            let mis2 = semi[i2].4.clone();
            let mut mis3 = Vec::<i32>::new();
            for l in l1 + len1..l2 {
                if b_seq[l as usize] != refs[t1 as usize].get((l + off1) as usize) {
                    mis3.push(l);
                }
            }
            let nmis = mis1.len() + mis2.len() + mis3.len();
            if nmis as f64 / ((l2 + len2) - l1) as f64 > MAX_RATE {
                continue;
            }
            semi[i1].3 = (l2 + len2) - l1;
            semi[i1].4.append(&mut mis3);
            semi[i1].4.append(&mut mis2.clone());
            unique_sort(&mut semi[i1].4);
            semi[i2].0 = -1_i32;
            to_delete[i2] = true;
        }
    }
    erase_if(&mut semi, &to_delete);
    report_semis(
        verbose,
        "SEMI ALIGNMENTS AFTER EXTENSION",
        &semi,
        &b_seq,
        refs,
        log,
    );

    // Merge overlapping alignments.
    // semi = {(t, off, pos on b, len, positions on b of mismatches)}

    let mut to_delete = vec![false; semi.len()];
    let mut i = 0;
    while i < semi.len() {
        let mut j = i + 1;
        while j < semi.len() {
            if semi[j].0 != semi[i].0 || semi[j].1 != semi[i].1 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in k1 + 1..j {
                if to_delete[k1] || to_delete[k2] {
                    continue;
                }
                let start1 = semi[k1].2;
                let start2 = semi[k2].2;
                let len1 = semi[k1].3;
                let len2 = semi[k2].3;
                let stop1 = start1 + len1;
                let stop2 = start2 + len2;
                let start = min(start1, start2);
                let stop = max(stop1, stop2);
                if stop - start <= len1 + len2 {
                    semi[k1].2 = start;
                    semi[k1].3 = stop - start;
                    let mut m2 = semi[k2].4.clone();
                    semi[k1].4.append(&mut m2);
                    unique_sort(&mut semi[k1].4);
                    to_delete[k2] = true;
                }
            }
        }
        i = j;
    }
    erase_if(&mut semi, &to_delete);
    report_semis(
        verbose,
        "SEMI ALIGNMENTS AFTER MERGER",
        &semi,
        &b_seq,
        refs,
        log,
    );

    // If a V gene aligns starting at 0, and goes at least 60% of the way to the end, and there
    // is only one alignment of the V gene, extend it to the end.
    // (Only one requirement ameliorated.)

    let mut i = 0;
    while i < semi.len() {
        let mut j = i + 1;
        while j < semi.len() {
            if semi[j].0 != semi[i].0 {
                break;
            }
            j += 1;
        }
        let mut k = i;
        let mut ok = false;
        if j - i == 1 {
            ok = true;
        } else if j - i == 2 {
            ok = true;
            if semi[i].2 < semi[i + 1].2 {
                k = i + 1;
            }
        }
        if ok {
            let offset = semi[k].1;
            let ref_start = semi[k].1 + semi[k].2;
            let tig_start = semi[k].2;
            let t = semi[k].0 as usize;
            if !rheaders[t].contains("segment") && refdata.is_v(t) {
                let r = &refs[t];
                let len = semi[k].3;
                if ref_start + len < r.len() as i32
                    && (ref_start + len) as f64 / r.len() as f64 >= 0.60
                    && len + tig_start < b_seq.len() as i32
                {
                    let start = ref_start + len;
                    let stop = min(r.len() as i32, b_seq.len() as i32 + offset);
                    for m in start..stop {
                        if b_seq[(m - offset) as usize] != r.get(m as usize) {
                            semi[k].4.push(m - offset);
                        }
                    }
                    semi[k].3 += stop - start;
                }
            }
        }
        i = j;
    }

    // Make sure that mismatches are unique sorted.

    for i in 0..semi.len() {
        unique_sort(&mut semi[i].4);
    }
    report_semis(
        verbose,
        "SEMI ALIGNMENTS AFTER SECOND EXTENSION",
        &semi,
        &b_seq,
        refs,
        log,
    );

    // Delete some subsumed alignments.

    let mut to_delete = vec![false; semi.len()];
    let mut i = 0;
    while i < semi.len() {
        let mut j = i + 1;
        while j < semi.len() {
            if semi[j].0 != semi[i].0 || semi[j].1 != semi[i].1 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in i..j {
                if semi[k1].1 + semi[k1].2 + semi[k1].3 == semi[k2].1 + semi[k2].2 + semi[k2].3
                    && semi[k1].3 > semi[k2].3
                {
                    to_delete[k2] = true;
                }
            }
        }
        i = j;
    }
    erase_if(&mut semi, &to_delete);
    report_semis(
        verbose,
        "SEMI ALIGNMENTS AFTER SUBSUMPTION",
        &semi,
        &b_seq,
        refs,
        log,
    );

    // Transform to create annx, having structure:
    // { ( sequence start, match length, ref tig, ref tig start, {mismatches} ) }.

    let mut annx = Vec::<(i32, i32, i32, i32, Vec<i32>)>::new();
    for x in semi.iter() {
        annx.push((x.2, x.3, x.0, x.2 + x.1, x.4.clone()));
    }
    unique_sort(&mut annx);

    // Delete matches that are 'too improper'.

    if !allow_improper {
        let mut to_delete: Vec<bool> = vec![false; annx.len()];
        for annxi in annx.iter_mut() {
            std::mem::swap(&mut annxi.0, &mut annxi.2);
            std::mem::swap(&mut annxi.1, &mut annxi.3);
        }
        annx.sort();
        let mut i1 = 0;
        loop {
            if i1 == annx.len() {
                break;
            }
            let j1 = next_diff1_5(&annx, i1 as i32);
            let mut min_imp = 1000000000;
            for k in i1..j1 as usize {
                let imp = min(annx[k].1, annx[k].2);
                min_imp = min(imp, min_imp);
            }
            const MAX_IMP: i32 = 60;
            if min_imp > MAX_IMP {
                for k in i1..j1 as usize {
                    to_delete[k] = true;
                }
            }
            i1 = j1 as usize;
        }
        erase_if(&mut annx, &to_delete);
        for annxi in annx.iter_mut() {
            std::mem::swap(&mut annxi.0, &mut annxi.2);
            std::mem::swap(&mut annxi.1, &mut annxi.3);
        }
        annx.sort();
    }

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nINITIAL ALIGNMENTS\n");
        for annxi in &annx {
            print_alignx(log, annxi, refdata);
        }
    }

    // Amongst V segments starting at zero on the V segment, if some start with
    // a start codon, delete the others.

    let mut have_starter = false;
    for annxi in &annx {
        let t = annxi.2 as usize;
        if !rheaders[t].contains("segment") && refdata.is_v(t) && annxi.3 == 0 {
            let p = annxi.0 as usize;
            if b_seq[p] == 0 // A
                && b_seq[p+1] == 3 // T
                && b_seq[p+2] == 2
            {
                // G
                have_starter = true;
                break;
            }
        }
    }
    if have_starter {
        let to_delete: Vec<bool> = annx
            .iter()
            .map(|annxi| {
                let t = annxi.2 as usize;
                if !rheaders[t].contains("segment") && refdata.is_v(t) && annxi.3 == 0 {
                    let p = annxi.0 as usize;
                    if !(b_seq[p] == 0 && b_seq[p + 1] == 3 && b_seq[p + 2] == 2) {
                        return true;
                    }
                }
                false
            })
            .collect();
        erase_if(&mut annx, &to_delete);
    }

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nALIGNMENTS ONE\n");
        for i in 0..annx.len() {
            print_alignx(log, &annx[i], refdata);
        }
    }

    // Remove inferior matches of the edge.  Two alignments are compared if the
    // length of their overlap on the contig is at least 85% of one of the alignment
    // lengths len1 and len2.  We compute the mismatch rates r1 and r2 between the
    // overlap interval and the respective references.  The first alignment wins if
    // at least one of the following happens:
    // 1. len1  > len2 and r1 <= r2
    // 2. len1 >= len2 and r2 <  r2
    // 3. len1 >= 1.5 * len2.
    //
    // Modified: multiple aligns of the same V segment are now group together in
    // the calculation.   And add indel penalty.
    //
    // ◼ For efficiency, inner loop should check to see if already deleted.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    let mut ts: Vec<(usize, usize)> = annx
        .iter()
        .enumerate()
        .map(|(i, annxi)| (annxi.2 as usize, i))
        .collect(); // { ( contig index, annx index ) }
    ts.sort_unstable();
    let mut i1 = 0;
    while i1 < ts.len() {
        let j1 = next_diff1_2(&ts, i1 as i32) as usize;
        let mut tlen1 = 0;
        for k in i1..j1 {
            tlen1 += annx[ts[k].1].1;
        }
        let mut i2 = 0;
        while i2 < ts.len() {
            let j2 = next_diff1_2(&ts, i2 as i32) as usize;
            let mut tlen2 = 0;
            for k in i2..j2 {
                tlen2 += annx[ts[k].1].1;
            }
            let (mut m1, mut m2) = (0, 0);
            let mut over = 0_i64;
            let mut offsets1 = Vec::<i32>::new();
            let mut offsets2 = Vec::<i32>::new();
            for k1 in i1..j1 {
                let u1 = ts[k1].1;
                offsets1.push(annx[u1].0 - annx[u1].3);
            }
            for k2 in i2..j2 {
                let u2 = ts[k2].1;
                offsets2.push(annx[u2].0 - annx[u2].3);
            }
            offsets1.sort_unstable();
            offsets2.sort_unstable();
            m1 += offsets1[offsets1.len() - 1] - offsets1[0];
            m2 += offsets2[offsets2.len() - 1] - offsets2[0];
            for k1 in i1..j1 {
                let u1 = ts[k1].1;
                let l1 = annx[u1].0;
                let len1 = annx[u1].1;
                for k2 in i2..j2 {
                    let u2 = ts[k2].1;
                    let l2 = annx[u2].0;
                    let len2 = annx[u2].1;
                    let start = max(l1, l2);
                    let stop = min(l1 + len1, l2 + len2);
                    if start >= stop {
                        continue;
                    }
                    over += stop as i64;
                    over -= start as i64;
                    for x in annx[u1].4.iter() {
                        if *x >= start && *x < stop {
                            m1 += 1;
                        }
                    }
                    for x in annx[u2].4.iter() {
                        if *x >= start && *x < stop {
                            m2 += 1;
                        }
                    }
                }
            }

            // Get mismatch rates.

            let (r1, r2) = (m1 as f64 / tlen1 as f64, m2 as f64 / tlen2 as f64);

            // Require that one of the intervals is at least 85% overlapped.

            const MIN_OVERLAP_FRAC: f64 = 0.85;
            if over as f64 / (min(tlen1, tlen2) as f64) >= MIN_OVERLAP_FRAC {
                // Decide if the second match is inferior.

                if (tlen1 > tlen2 && r1 <= r2)
                    || (tlen1 >= tlen2 && r1 < r2)
                    || tlen1 as f64 >= 1.5 * tlen2 as f64
                {
                    if verbose {
                        fwriteln!(
                            log,
                            "\nsee tlen1 = {}, tlen2 = {}, m1 = {}, m2 = {}, \
                             r1 = {:.3}, r2 = {:.3}\nthis alignment",
                            tlen1,
                            tlen2,
                            m1,
                            m2,
                            r1,
                            r2
                        );
                        for k in i1..j1 {
                            let t = ts[k].0;
                            let u = ts[k].1;
                            let l = annx[u].0;
                            let len = annx[u].1;
                            let p = annx[u].3;
                            let mis = annx[u].4.len();
                            fwriteln!(
                                log,
                                "{}-{} ==> {}-{} on {}(mis={})",
                                l,
                                l + len,
                                p,
                                p + len,
                                rheaders[t],
                                mis
                            );
                        }
                        fwriteln!(log, "beats this alignment");
                        for k in i2..j2 {
                            let t = ts[k].0;
                            let u = ts[k].1;
                            let l = annx[u].0;
                            let len = annx[u].1;
                            let p = annx[u].3;
                            let mis = annx[u].4.len();
                            fwriteln!(
                                log,
                                "{}-{} ==> {}-{} on {}(mis={})",
                                l,
                                l + len,
                                p,
                                p + len,
                                rheaders[t],
                                mis
                            );
                        }
                    }
                    for k in i2..j2 {
                        to_delete[ts[k].1] = true;
                    }
                }
            }
            i2 = j2;
        }
        i1 = j1;
    }
    erase_if(&mut annx, &to_delete);

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nALIGNMENTS TWO\n");
        for i in 0..annx.len() {
            print_alignx(log, &annx[i], refdata);
        }
    }

    // If there are two alignments to a particular V region, or a UTR, try to edit
    // them so that their start/stop position abut perfect on one side (either the
    // contig or the reference), and do not overlap on the other side, thus
    // exhibiting an indel.
    // ◼ The approach to answering this seems very inefficient.
    // ◼ When this was moved here, some UTR alignments disappeared.

    // { ( sequence start, match length, ref tig, ref tig start, {mismatches} ) }.
    if abut {
        let mut to_delete: Vec<bool> = vec![false; annx.len()];
        let mut aligns = vec![0; refs.len()];
        for i in 0..annx.len() {
            aligns[annx[i].2 as usize] += 1;
        }
        for i1 in 0..annx.len() {
            let t = annx[i1].2 as usize;
            if rheaders[t].contains("segment") {
                continue;
            }
            if !rheaders[t].contains("segment") && !refdata.is_u(t) && !refdata.is_v(t) {
                continue;
            }
            let off1 = annx[i1].3 - annx[i1].0;
            for i2 in 0..annx.len() {
                if to_delete[i1] || to_delete[i2] {
                    continue;
                }
                if i2 == i1 || annx[i2].2 as usize != t {
                    continue;
                }
                let (l1, mut l2) = (annx[i1].0 as usize, annx[i2].0 as usize);
                if l1 >= l2 {
                    continue;
                }
                let (mut len1, mut len2) = (annx[i1].1 as usize, annx[i2].1 as usize);
                if l1 + len1 > l2 + len2 {
                    continue;
                }
                let (p1, p2) = (annx[i1].3, annx[i2].3);
                let (start1, stop1) = (l1, (l2 + len2)); // extent on contig
                let (start2, stop2) = (p1 as usize, (p2 as usize + len2)); // extent on ref
                if !(start1 < stop1 && start2 < stop2) {
                    continue;
                }
                let tot1 = stop1 - start1;
                let tot2 = stop2 - start2;
                let off2 = annx[i2].3 - annx[i2].0;

                // Case where there is no indel.

                if tot1 == tot2 && aligns[t] == 2 {
                    let mut mis = annx[i1].4.clone();
                    for p in l1 + len1..l2 {
                        if b_seq[p] != refs[t].get((p as i32 + off1) as usize) {
                            mis.push(p as i32);
                        }
                    }
                    mis.append(&mut annx[i2].4.clone());
                    annx[i1].1 = tot1 as i32;
                    annx[i1].4 = mis;
                    to_delete[i2] = true;
                    continue;
                }

                // Case of insertion.

                if tot1 > tot2 && aligns[t] == 2 {
                    let start1 = start1 as i32;
                    let stop1 = stop1 as i32;
                    let ins = (tot1 - tot2) as i32;
                    let mut best_ipos = 0;
                    let mut best_mis = 1000000;
                    let mut best_mis1 = Vec::<i32>::new();
                    let mut best_mis2 = Vec::<i32>::new();
                    for ipos in start1..=stop1 - ins {
                        let mut mis1 = Vec::<i32>::new();
                        let mut mis2 = Vec::<i32>::new();
                        for p in start1..ipos {
                            if b_seq[p as usize] != refs[t].get((p + off1) as usize) {
                                mis1.push(p);
                            }
                        }
                        for p in ipos + ins..stop1 {
                            if b_seq[p as usize] != refs[t].get((p + off2) as usize) {
                                mis2.push(p);
                            }
                        }
                        let mis = (mis1.len() + mis2.len()) as i32;
                        if mis < best_mis {
                            best_mis = mis;
                            best_mis1 = mis1;
                            best_mis2 = mis2;
                            best_ipos = ipos;
                        }
                    }
                    annx[i1].1 = best_ipos - start1;
                    annx[i1].4 = best_mis1;
                    annx[i2].1 = stop1 - best_ipos - ins;
                    annx[i2].0 = best_ipos + ins;
                    annx[i2].3 = best_ipos + ins + off2;
                    annx[i2].4 = best_mis2;
                    continue;
                }

                // Case of deletion.

                if tot1 < tot2 && aligns[t] == 2 {
                    let start2 = start2 as i32;
                    let stop2 = stop2 as i32;
                    let del = (tot2 - tot1) as i32;
                    let mut best_dpos = 0;
                    let mut best_mis = 1000000;
                    let mut best_mis1 = Vec::<i32>::new();
                    let mut best_mis2 = Vec::<i32>::new();
                    for dpos in start2..=stop2 - del {
                        let mut mis1 = Vec::<i32>::new();
                        let mut mis2 = Vec::<i32>::new();
                        for q in start2..dpos {
                            let p = q - off1;
                            if b_seq[p as usize] != refs[t].get(q as usize) {
                                mis1.push(p);
                            }
                        }
                        for q in dpos + del..stop2 {
                            let p = q - off2;
                            if b_seq[p as usize] != refs[t].get(q as usize) {
                                mis2.push(p);
                            }
                        }
                        let mis = (mis1.len() + mis2.len()) as i32;
                        if mis < best_mis {
                            best_mis = mis;
                            best_mis1 = mis1;
                            best_mis2 = mis2;
                            best_dpos = dpos;
                        }
                    }
                    annx[i1].1 = best_dpos - start2;
                    annx[i1].4 = best_mis1;
                    annx[i2].0 = best_dpos + del - off2;
                    annx[i2].1 = stop2 - best_dpos - del;
                    annx[i2].3 = best_dpos + del;
                    annx[i2].4 = best_mis2;
                    continue;
                }

                // It's not clear why the rest of this code helps, but it does.

                let p1 = p1 as usize;
                let mut p2 = p2 as usize;

                let b1 = b.slice(start1, stop1).to_owned();
                let b2 = refs[t].slice(start2, stop2).to_owned();

                let a = affine_align(&b1, &b2);
                let mut del = Vec::<(usize, usize, usize)>::new();
                let mut ins = Vec::<(usize, usize, usize)>::new();
                let ops = &a.operations;
                let mut i = 0;
                let (mut z1, mut z2) = (l1 + a.xstart, p1 + a.ystart);
                if a.ystart > 0 {
                    continue;
                }
                let mut matches = 0;
                while i < ops.len() {
                    let mut opcount = 1;
                    while i + opcount < ops.len()
                        && (ops[i] == Del || ops[i] == Ins)
                        && ops[i] == ops[i + opcount]
                    {
                        opcount += 1;
                    }
                    match ops[i] {
                        Match => {
                            matches += 1;
                            z1 += 1;
                            z2 += 1;
                        }
                        Subst => {
                            z1 += 1;
                            z2 += 1;
                        }
                        Del => {
                            del.push((z1, z2, opcount));
                            if verbose {
                                fwriteln!(log, "\nsee del[{}]", opcount);
                            }
                            z2 += opcount;
                        }
                        Ins => {
                            ins.push((z1, z2, opcount));
                            if verbose {
                                fwriteln!(log, "\nsee ins[{}]", opcount);
                            }
                            z1 += opcount;
                        }
                        Xclip(d) => {
                            z1 += d;
                        }
                        Yclip(d) => {
                            z2 += d;
                        }
                    }
                    i += opcount;
                }
                if verbose {
                    fwriteln!(log, "\ntrying to merge\n{}\n{}", rheaders[t], rheaders[t]);
                    fwriteln!(log, "|del| = {}, |ins| = {}", del.len(), ins.len());
                }
                if del.solo() && ins.is_empty() {
                    let (l, p, n) = (del[0].0, del[0].1, del[0].2);
                    if n != (p2 + len2 - p1) - (l2 + len2 - l1) {
                        continue;
                    }
                    len1 = l - l1;
                    if len1 >= matches {
                        continue;
                    }
                    len2 = l2 + len2 - l1 - len1;
                    l2 = l;
                    p2 = p + n;
                }
                if del.is_empty() && ins.solo() {
                    let (l, p, n) = (ins[0].0, ins[0].1, ins[0].2);
                    // ◼ This is buggy.  It fails if overflow detection is on.
                    if n + l1 + p2 != l2 + p1 {
                        continue;
                    }
                    len1 = l - l1;
                    if len1 >= matches {
                        continue;
                    }
                    len2 = p2 + len2 - p1 - len1;
                    l2 = l + n;
                    p2 = p;
                }
                if del.len() + ins.len() == 0 {
                    to_delete[i2] = true;
                    len1 = (annx[i2].0 + annx[i2].1 - annx[i1].0) as usize;
                    annx[i1].1 = len1 as i32;
                    annx[i1].4.truncate(0);
                    for j in 0..len1 {
                        if b_seq[l1 + j] != refs[t].get(p1 + j) {
                            annx[i1].4.push((l1 + j) as i32);
                        }
                    }
                }
                if del.len() + ins.len() == 1 {
                    annx[i2].0 = l2 as i32;
                    annx[i1].1 = len1 as i32;
                    annx[i2].1 = len2 as i32;
                    annx[i2].3 = p2 as i32;
                    annx[i1].4.truncate(0);
                    annx[i2].4.truncate(0);
                    for j in 0..len1 {
                        if b_seq[l1 + j] != refs[t].get(p1 + j) {
                            annx[i1].4.push((l1 + j) as i32);
                        }
                    }
                    for j in 0..len2 {
                        if b_seq[l2 + j] != refs[t].get(p2 + j) {
                            annx[i2].4.push((l2 + j) as i32);
                        }
                    }
                }
            }
        }
        erase_if(&mut annx, &to_delete);
    }

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nALIGNMENTS THREE\n");
        for i in 0..annx.len() {
            print_alignx(log, &annx[i], refdata);
        }
    }

    // Choose between segments if one clearly wins.  For this calculation, we
    // put UTR and V segments together.  The way the choice is made could be refined.
    //
    // ◼ At least in some cases, a better way of comparing errors would be to first
    // ◼ extend the alignments so that their endpoints on the contigs agree, to the
    // ◼ extent that this is possible.
    //
    // ◼ This code has not really been adequately tested to see if the right
    // ◼ choices are being made.
    //
    // ◼ Note danger with nonstandard references.
    //
    // ◼ Really should have ho_interval here.

    let mut combo = Vec::<(String, i32, usize)>::new();
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if !rheaders[t].contains("segment") {
            combo.push((
                refdata.name[t].clone() + "." + &refdata.transcript[t],
                refdata.id[t],
                i,
            ));
        }
    }
    combo.sort();
    //                     cov                 mis    locs        rstarts     mis_nutr
    let mut data = Vec::<(Vec<(usize, usize)>, usize, Vec<usize>, Vec<usize>, usize)>::new();
    let mut i = 0;
    while i < combo.len() {
        let j = next_diff1_3(&combo, i as i32) as usize;
        let mut cov = Vec::<(usize, usize)>::new();
        let mut mis = 0;
        let mut mis_nutr = 0;
        let mut locs = Vec::<usize>::new();
        let mut rstarts = Vec::<usize>::new();
        for k in i..j {
            locs.push(combo[k].2);
            let a = &annx[combo[k].2];
            rstarts.push(a.3 as usize);
            cov.push((a.0 as usize, (a.0 + a.1) as usize));
            mis += a.4.len();
            if !refdata.is_u(a.2 as usize) {
                mis_nutr += a.4.len();
            }
        }
        data.push((cov, mis, locs, rstarts, mis_nutr));
        i = j;
    }
    let mut to_delete = vec![false; annx.len()];
    let mut deleted = vec![false; data.len()];
    for i1 in 0..data.len() {
        if deleted[i1] {
            continue;
        }
        for i2 in 0..data.len() {
            if i2 == i1 {
                continue;
            }
            let t1 = annx[data[i1].2[0]].2 as usize;
            let t2 = annx[data[i2].2[0]].2 as usize;

            let same_class = (refdata.segtype[t1] == refdata.segtype[t2])
                || (refdata.is_v(t1) && refdata.is_u(t2))
                || (refdata.is_u(t1) && refdata.is_v(t2));
            if !same_class {
                continue;
            }

            // At this point we have two alignments, for which the reference sequence names
            // are the same, and for which either they both have the same segment type
            // (U, V, J, C), or else one is U and one is V.  Next we compare them, first
            // gathering information about UTRs.

            let name1 = &refdata.name[t1];
            let name2 = &refdata.name[t2];
            let (mut utr1, mut utr2) = (false, false);
            if refdata.is_v(t1) || refdata.is_u(t1) {
                utr1 = refdata.has_utr[name1];
                utr2 = refdata.has_utr[name2];
            }
            let (mut have_utr_align1, mut have_utr_align2) = (false, false);
            for j in data[i1].2.iter() {
                if refdata.is_u(annx[*j].2 as usize) {
                    have_utr_align1 = true;
                }
            }
            for j in data[i2].2.iter() {
                if refdata.is_u(annx[*j].2 as usize) {
                    have_utr_align2 = true;
                }
            }

            // Now find their mismatch positions.

            let n = b_seq.len();
            let (mut mis1, mut mis2) = (vec![false; n], vec![false; n]);
            for j in data[i1].2.iter() {
                for p in annx[*j].4.iter() {
                    mis1[*p as usize] = true;
                }
            }
            for j in data[i2].2.iter() {
                for p in annx[*j].4.iter() {
                    mis2[*p as usize] = true;
                }
            }

            // Compute the fraction of i2 coverage that's outside i1 coverage.
            // ◼ This is horrendously inefficient.  Use ho intervals.

            let (mut cov1, mut cov2) = (vec![false; n], vec![false; n]);
            for j in 0..data[i1].0.len() {
                let t = annx[data[i1].2[j]].2;
                if utr2 || !refdata.is_u(t as usize) {
                    let x = &data[i1].0[j];
                    for m in x.0..x.1 {
                        cov1[m] = true;
                    }
                }
            }
            for j in 0..data[i2].0.len() {
                let t = annx[data[i2].2[j]].2;
                if utr1 || !refdata.is_u(t as usize) {
                    let x = &data[i2].0[j];
                    for m in x.0..x.1 {
                        cov2[m] = true;
                    }
                }
            }
            let (mut total1, mut total2) = (0, 0);
            for l in 0..n {
                if cov1[l] {
                    total1 += 1;
                }
            }
            for l in 0..n {
                if cov2[l] {
                    total2 += 1;
                }
            }

            // Same as above but always exclude UTRs.

            let (mut cov1_nu, mut cov2_nu) = (vec![false; n], vec![false; n]);
            for j in 0..data[i1].0.len() {
                let t = annx[data[i1].2[j]].2;
                if !refdata.is_u(t as usize) {
                    let x = &data[i1].0[j];
                    for m in x.0..x.1 {
                        cov1_nu[m] = true;
                    }
                }
            }
            for j in 0..data[i2].0.len() {
                let t = annx[data[i2].2[j]].2;
                if !refdata.is_u(t as usize) {
                    let x = &data[i2].0[j];
                    for m in x.0..x.1 {
                        cov2_nu[m] = true;
                    }
                }
            }
            let (mut total1_nu, mut total2_nu) = (0, 0);
            for l in 0..n {
                if cov1_nu[l] {
                    total1_nu += 1;
                }
            }
            for l in 0..n {
                if cov2_nu[l] {
                    total2_nu += 1;
                }
            }

            // Compute amount shared.

            let mut share = 0;
            for l in 0..n {
                if cov1[l] && cov2[l] {
                    share += 1;
                }
            }
            let outside1 = percent_ratio(total1 - share, total1);
            let outside2 = percent_ratio(total2 - share, total2);

            // Find the number of mismatches in the overlap region.

            let (mut m1, mut m2) = (0, 0);
            for l in 0..n {
                if cov1[l] && cov2[l] {
                    if mis1[l] {
                        m1 += 1;
                    }
                    if mis2[l] {
                        m2 += 1;
                    }
                }
            }

            // Compute error rates.
            // ◼ This is incorrect in the case where the UTR has been excluded.
            // Added separate estimates with UTRs excluded.

            let err1 = percent_ratio(data[i1].1, total1);
            let err2 = percent_ratio(data[i2].1, total2);
            let err1_nu = percent_ratio(data[i1].4, total1_nu);
            let err2_nu = percent_ratio(data[i2].4, total2_nu);

            // Compute zstops.

            let (mut zstop1, mut zstop2) = (0, 0);
            for l in 0..data[i1].2.len() {
                let t = annx[data[i1].2[l]].2 as usize;
                if refdata.is_v(t) && (data[i1].3[l] == 0 || data[i1].0[l].0 == 0) {
                    zstop1 = max(zstop1, data[i1].0[l].1);
                }
            }
            for l in 0..data[i2].2.len() {
                let t = annx[data[i2].2[l]].2 as usize;
                if refdata.is_v(t) && (data[i2].3[l] == 0 || data[i2].0[l].0 == 0) {
                    zstop2 = max(zstop2, data[i2].0[l].1);
                }
            }

            // Decide if the first wins.  And symmetrize to prevent double
            // deletion.  Be very careful to respect this if editing!

            let (mut win1, mut win2) = (false, false);
            let c1 = m1 == m2 && !have_utr_align2 && err1_nu == err2_nu && outside1 > outside2;
            let c2 = m2 == m1 && !have_utr_align1 && err2_nu == err1_nu && outside2 > outside1;
            if (zstop1 > zstop2 + 20 && (outside2 <= 10.0 || total2 - share <= 10))
                || (outside1 >= 10.0 && outside2 <= 1.0 && err1 - err2 <= 2.5)
            {
                win1 = true;
            } else if zstop1 == 0 && zstop2 > 0 {
            } else if (outside2 <= 10.0 || total2 - share <= 10)
                && (m1 < m2
                    || (m1 == m2 && err1 < err2 && !c2)
                    || (m1 == m2 && err1 == err2 && outside1 > outside2)
                    || (m1 == m2 && err1 == err2 && outside1 == outside2 && t1 < t2)
                    || c1)
            {
                win1 = true;
            }

            // Symmetrization.

            if (zstop2 > zstop1 + 20 && (outside1 <= 10.0 || total1 - share <= 10))
                || (outside2 >= 10.0 && outside1 <= 1.0 && err2 - err1 <= 2.5)
            {
                win2 = true;
            } else if zstop2 == 0 && zstop1 > 0 {
            } else if (outside1 <= 10.0 || total1 - share <= 10)
                && (m2 < m1
                    || (m2 == m1 && err2 < err1 && !c1)
                    || (m2 == m1 && err2 == err1 && outside2 > outside1)
                    || (m2 == m1 && err2 == err1 && outside2 == outside1 && t2 < t1)
                    || c2)
            {
                win2 = true;
            }
            if win2 {
                win1 = false;
            }

            // Verbose logging.

            if verbose {
                fwriteln!(log, "\nCOMPARING");
                for l in 0..data[i1].2.len() {
                    let t = annx[data[i1].2[l]].2 as usize;
                    let cov = &data[i1].0[l];
                    let mis = data[i1].1;
                    fwriteln!(
                        log,
                        "{}, cov = {}-{}, mis = {}",
                        rheaders[t],
                        cov.0,
                        cov.1,
                        mis
                    );
                }
                fwriteln!(log, "TO");
                for l in 0..data[i2].2.len() {
                    let t = annx[data[i2].2[l]].2 as usize;
                    let cov = &data[i2].0[l];
                    let mis = data[i2].1;
                    fwriteln!(
                        log,
                        "{}, cov = {}-{}, mis = {}",
                        rheaders[t],
                        cov.0,
                        cov.1,
                        mis
                    );
                }
                fwriteln!(log, "zstop1 = {}, zstop2 = {}", zstop1, zstop2);
                fwriteln!(log, "m1 = {}, m2 = {}", m1, m2);
                fwriteln!(log, "err1 = {}, err2 = {}", err1, err2);
                fwriteln!(log, "err1_nu = {}, err2_nu = {}", err1_nu, err2_nu);
                fwriteln!(log, "utr1 = {}, utr2 = {}", utr1, utr2);
                fwriteln!(
                    log,
                    "have_utr_align1 = {}, have_utr_align2 = {}",
                    have_utr_align1,
                    have_utr_align2
                );
                fwriteln!(
                    log,
                    "total1 = {}, total2 = {}, share = {}",
                    total1,
                    total2,
                    share
                );
                fwriteln!(
                    log,
                    "outside1 = {:.1}%, outside2 = {:.1}%, \
                     total2 - share = {}, err1 = {:.1}%, err2 = {:.1}%",
                    outside1,
                    outside2,
                    total2 - share,
                    err1,
                    err2
                );
                fwriteln!(log, "win1 = {}, win2 = {}", win1, win2);
            }

            // Pick "randomly" in case of tie.  Unless we're comparing TRBC1 to TRBC2, and in
            // that case, do nothing.  This is needed because of the handling below of these
            // segments.

            if outside1 == 0.0
                && outside2 == 0.0
                && zstop1 == zstop2
                && m1 == m2
                && err1 == err2
                && t1 < t2
            {
                if refdata.name[t1] == *"TRBC1" && refdata.name[t2] == *"TRBC2" {
                    continue;
                }
                if refdata.name[t2] == *"TRBC1" && refdata.name[t1] == *"TRBC2" {
                    continue;
                }
                win1 = true;
            }

            // Make decision.

            if win1 {
                for l in data[i2].2.iter() {
                    to_delete[*l] = true;
                }
                deleted[i2] = true;
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Delete some subsumed alignments.

    let mut to_delete = vec![false; annx.len()];
    let mut i = 0;
    while i < annx.len() {
        let mut j = i + 1;
        while j < annx.len() {
            if annx[j].0 != annx[i].0 {
                break;
            }
            j += 1;
        }
        for k1 in i..j {
            for k2 in i..j {
                if annx[k1].2 == annx[k2].2 && annx[k1].3 == annx[k2].3 && annx[k1].1 > annx[k2].1 {
                    to_delete[k2] = true;
                }
            }
        }
        i = j;
    }
    erase_if(&mut annx, &to_delete);

    // Extend some alignments.
    // { ( sequence start, match length, ref tig, ref tig start, {mismatches} ) }.

    let mut aligns = vec![0; refs.len()];
    for i in 0..annx.len() {
        aligns[annx[i].2 as usize] += 1;
    }
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        let len = annx[i].1 as usize;
        if aligns[t] == 1
            && annx[i].3 == 0
            && len < refs[t].len()
            && len as f64 / refs[t].len() as f64 >= 0.75
            && (refs[t].len() as i32 + annx[i].0 - annx[i].3) as usize <= b_seq.len()
        {
            for p in len..refs[t].len() {
                let q = p as i32 + annx[i].0 - annx[i].3;
                if b_seq[q as usize] != refs[t].get(p) {
                    annx[i].4.push(q);
                }
            }
            annx[i].1 = refs[t].len() as i32;
        }
    }

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nALIGNMENTS FOUR\n");
        for i in 0..annx.len() {
            print_alignx(log, &annx[i], refdata);
        }
    }

    // If two V segments are aligned starting at 0 on the reference and one
    // is aligned a lot further, it wins.

    let mut lens = vec![0; refdata.refs.len()];
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        lens[t] += annx[i].3 + annx[i].1;
    }
    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        for i2 in 0..annx.len() {
            let (t1, t2) = (annx[i1].2 as usize, annx[i2].2 as usize);
            if rheaders[t1].contains("segment") || rheaders[t2].contains("segment") {
                continue;
            }
            if !refdata.is_v(t1) || !refdata.is_v(t2) {
                continue;
            }
            if t1 == t2 {
                continue;
            }
            let (p1, p2) = (annx[i1].3, annx[i2].3);
            if p1 > 0 {
                continue;
            }
            const MIN_EXT: i32 = 50;
            if (p2 > 0 && lens[t1] >= lens[t2]) || (p2 == 0 && lens[t1] >= lens[t2] + MIN_EXT) {
                if verbose {
                    fwriteln!(log, "");
                    print_alignx(log, &annx[i1], refdata);
                    fwriteln!(log, "beats");
                    print_alignx(log, &annx[i2], refdata);
                }
                to_delete[i2] = true;
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // For IG, if we have a C segment that aligns starting at zero, and a V segment
    // that aligns, but no J segment, try to find a J segment alignment.  For now we
    // assume that the J aligns up to exactly the point where the C starts, or to
    // one base after.  We require that the last 20 bases of the J match with at
    // most 5 mismatches.

    let (mut igv, mut igj) = (false, false);
    let mut igc = -1_i32;
    const J_TOT: i32 = 20;
    const J_MIS: i32 = 5;
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if rheaders[t].contains("segment") {
            continue;
        }
        let rt = refdata.rtype[t];
        if (0..3).contains(&rt) {
            if refdata.segtype[t] == "V" {
                igv = true;
            } else if refdata.segtype[t] == "J" {
                igj = true;
            } else if refdata.segtype[t] == "C"
                && annx[i].3 == 0
                && annx[i].0 >= J_TOT
                && refs[t].len() >= J_TOT as usize
            {
                igc = annx[i].0;
            }
        }
    }
    if igc >= 0 && igv && !igj {
        let mut best_t = -1_i32;
        let mut best_mis = 1000000;
        let mut best_z = -1_i32;
        for z in 0..2 {
            for l in 0..refdata.igjs.len() {
                let t = refdata.igjs[l];
                let n = refs[t].len();
                if n > igc as usize + z {
                    continue;
                }
                let i = igc as usize + z - n; // start of J on contig
                let (mut total, mut mis) = (0, 0);
                for j in (0..n).rev() {
                    total += 1;
                    if b_seq[i + j] != refs[t].get(j) {
                        mis += 1;
                        if total <= J_TOT && mis > J_MIS {
                            break;
                        }
                    }
                }
                if total == n as i32 && mis < best_mis {
                    best_t = t as i32;
                    best_mis = mis;
                    best_z = z as i32;
                }
            }
        }
        if best_t >= 0 {
            let t = best_t as usize;
            let n = refs[t].len() as i32;
            let i = igc + best_z - n;
            let mut mis = Vec::<i32>::new();
            for j in 0..n {
                if b_seq[(i + j) as usize] != refs[t].get(j as usize) {
                    mis.push(i + j);
                }
            }
            annx.push((i, n, best_t, 0_i32, mis));
            annx.sort();
        }
    }

    // If there is a D gene alignment that is from a different chain type than the V gene
    // alignment, delete it.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        let t1 = annx[i1].2 as usize;
        if !rheaders[t1].contains("segment") && refdata.segtype[t1] == "D" {
            let mut have_v = false;
            for i2 in 0..annx.len() {
                let t2 = annx[i2].2 as usize;
                if !rheaders[t2].contains("segment")
                    && refdata.segtype[t2] == "V"
                    && refdata.rtype[t1] == refdata.rtype[t2]
                {
                    have_v = true;
                }
            }
            if !have_v {
                to_delete[i1] = true;
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // For IGH and TRB (and TRD in Gamma/delta mode), if there is a V and J, but no D, look for a D that matches nearly perfectly
    // between them.  We consider only alignments having no indels.  The following conditions
    // are required:
    // 1. At most three mismatches.
    // 2. Excluding genes having the same name:
    //    (a) all others have more mismatches
    //    (b) all others have no more matches.

    let (mut v, mut d, mut j) = (false, false, false);
    let (mut vstop, mut jstart) = (0, 0);
    const VJTRIM: i32 = 10;
    let mut v_rtype = -2_i32;
    for ann in &annx {
        let t = ann.2 as usize;
        if !rheaders[t].contains("segment") {
            let rt = refdata.rtype[t];
            // IGH or TRB (or TRD in Gamma/delta mode)
            if rt == 0 || rt == 4 || rt == 5 {
                if refdata.segtype[t] == "V" {
                    v = true;
                    vstop = ann.0 + ann.1;
                    v_rtype = rt;
                } else if refdata.segtype[t] == "D" {
                    d = true;
                } else if refdata.segtype[t] == "J" {
                    j = true;
                    jstart = ann.0;
                }
            }
        }
    }
    if v && !d && j {
        let mut results = Vec::<(usize, usize, usize, String, usize, Vec<i32>)>::new();
        let start = max(0, vstop - VJTRIM);
        let stop = min(b.len() as i32, jstart + VJTRIM);
        const MAX_MISMATCHES: usize = 3;
        for t in refdata.ds.iter() {
            if refdata.rtype[*t] == v_rtype {
                let r = &refdata.refs[*t];
                for m in start..=stop - (r.len() as i32) {
                    let mut mismatches = Vec::<i32>::new();
                    for x in 0..r.len() {
                        if r.get(x) != b_seq[(m + x as i32) as usize] {
                            mismatches.push(x as i32);
                        }
                    }
                    let matches = r.len() - mismatches.len();
                    let mut gene = refdata.name[*t].clone();
                    if gene.contains('*') {
                        gene = gene.before("*").to_string();
                    }
                    results.push((mismatches.len(), matches, *t, gene, m as usize, mismatches));
                }
            }
        }
        results.sort();
        if !results.is_empty() && results[0].0 <= MAX_MISMATCHES {
            let mut to_delete = vec![false; results.len()];
            for i in 1..results.len() {
                if results[i].3 == results[0].3 {
                    to_delete[i] = true;
                }
            }
            erase_if(&mut results, &to_delete);
            if results.solo() || results[0].0 < results[1].0 {
                let mut best_matches = 0;
                for i in 0..results.len() {
                    best_matches = max(best_matches, results[i].1);
                }
                if results[0].1 == best_matches {
                    let t = results[0].2;
                    let r = results[0].0 + results[0].1;
                    let m = results[0].4;
                    annx.push((m as i32, r as i32, t as i32, 0, results[0].5.clone()));
                    annx.sort();
                }
            }
        }
    }

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nALIGNMENTS FIVE\n");
        for i in 0..annx.len() {
            print_alignx(log, &annx[i], refdata);
        }
    }

    // A J segment that goes up to its end beats any J segment that doesn't.
    // If they both go up to the end, choose.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        for i2 in 0..annx.len() {
            let (t1, t2) = (annx[i1].2 as usize, annx[i2].2 as usize);
            if rheaders[t1].contains("segment") || rheaders[t2].contains("segment") {
                continue;
            }
            if !refdata.is_j(t1) || !refdata.is_j(t2) {
                continue;
            }
            let (len1, len2) = (annx[i1].1, annx[i2].1);
            let (l1, l2) = (annx[i1].0, annx[i2].0);
            let (p1, p2) = (annx[i1].3, annx[i2].3);
            if len1 + p1 == refs[t1].len() as i32 && len2 + p2 < refs[t2].len() as i32 {
                to_delete[i2] = true;
            }
            if len1 + p1 == refs[t1].len() as i32 && len2 + p2 == refs[t2].len() as i32 {
                let (mut mis1, mut mis2) = (0, 0);
                let mut y1 = refs[t1].len() as i32 - 1;
                let mut y2 = refs[t2].len() as i32 - 1;
                let (mut x1, mut x2) = (y1 + l1 - p1, y2 + l2 - p2);
                loop {
                    if b_seq[x1 as usize] != refs[t1].get(y1 as usize) {
                        mis1 += 1;
                    }
                    if b_seq[x2 as usize] != refs[t2].get(y2 as usize) {
                        mis2 += 1;
                    }
                    if x1 == 0 || y1 == 0 || x2 == 0 || y2 == 0 {
                        break;
                    }
                    x1 -= 1;
                    y1 -= 1;
                    x2 -= 1;
                    y2 -= 1;
                }
                if mis1 < mis2 || (mis1 == mis2 && t1 < t2) {
                    to_delete[i2] = true;
                }
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Pick between C segments starting at zero.  And favor zero.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        for i2 in 0..annx.len() {
            if i2 == i1 {
                continue;
            }
            let (t1, t2) = (annx[i1].2 as usize, annx[i2].2 as usize);
            if rheaders[t1].contains("segment") || rheaders[t2].contains("segment") {
                continue;
            }
            if !refdata.is_c(t1) || !refdata.is_c(t2) {
                continue;
            }
            let (l1, l2) = (annx[i1].0, annx[i2].0);
            let (p1, p2) = (annx[i1].3, annx[i2].3);
            if p1 > 0 {
                continue;
            }
            if p1 == 0 && p2 > 0 {
                to_delete[i2] = true;
            }
            let (mut mis1, mut mis2) = (0, 0);
            let (mut y1, mut y2) = (p1, p2);
            let (mut x1, mut x2) = (l1, l2);
            loop {
                if b_seq[x1 as usize] != refs[t1].get(y1 as usize) {
                    mis1 += 1;
                }
                if b_seq[x2 as usize] != refs[t2].get(y2 as usize) {
                    mis2 += 1;
                }
                x1 += 1;
                y1 += 1;
                x2 += 1;
                y2 += 1;
                if x1 == b.len() as i32 || y1 == refs[t1].len() as i32 {
                    break;
                }
                if x2 == b.len() as i32 || y2 == refs[t2].len() as i32 {
                    break;
                }
            }

            if mis1 == mis2 {
                if refdata.name[t1] == *"TRBC1" && refdata.name[t2] == *"TRBC2" {
                    continue;
                }
                if refdata.name[t2] == *"TRBC1" && refdata.name[t1] == *"TRBC2" {
                    continue;
                }
            }
            if mis1 < mis2 || (mis1 == mis2 && t1 < t2) {
                to_delete[i2] = true;
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Pick between V segments starting at zero.  And favor zero.

    let mut nv = 0;
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if rheaders[t].contains("segment") {
            continue;
        }
        if refdata.is_v(t) {
            nv += 1;
        }
    }
    if nv == 2 {
        let mut to_delete: Vec<bool> = vec![false; annx.len()];
        for i1 in 0..annx.len() {
            for i2 in 0..annx.len() {
                let (t1, t2) = (annx[i1].2 as usize, annx[i2].2 as usize);
                if t2 == t1 {
                    continue;
                }
                if rheaders[t1].contains("segment") || rheaders[t2].contains("segment") {
                    continue;
                }
                if !refdata.is_v(t1) || !refdata.is_v(t2) {
                    continue;
                }
                let (l1, l2) = (annx[i1].0, annx[i2].0);
                let (p1, p2) = (annx[i1].3, annx[i2].3);
                if p1 > 0 {
                    continue;
                }
                if p2 > 0 {
                    to_delete[i2] = true;
                }
                let (mut mis1, mut mis2) = (0, 0);
                let (mut y1, mut y2) = (p1, p2);
                let (mut x1, mut x2) = (l1, l2);
                loop {
                    if b_seq[x1 as usize] != refs[t1].get(y1 as usize) {
                        mis1 += 1;
                    }
                    if b_seq[x2 as usize] != refs[t2].get(y2 as usize) {
                        mis2 += 1;
                    }
                    x1 += 1;
                    y1 += 1;
                    x2 += 1;
                    y2 += 1;
                    if x1 == b.len() as i32 || y1 == refs[t1].len() as i32 {
                        break;
                    }
                    if x2 == b.len() as i32 || y2 == refs[t2].len() as i32 {
                        break;
                    }
                }
                if mis1 < mis2 {
                    to_delete[i2] = true;
                }
            }
        }
        erase_if(&mut annx, &to_delete);
    }

    // Remove UTR annotations that have no matching V annotation.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    let (mut u, mut v) = (Vec::<String>::new(), Vec::<String>::new());
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if !rheaders[t].contains("segment") {
            let name = rheaders[t].after("|").between("|", "|");
            if rheaders[t].contains("UTR") {
                u.push(name.to_string());
            }
            if rheaders[t].contains("V-REGION") {
                v.push(name.to_string());
            }
        }
    }
    v.sort();
    for i in 0..u.len() {
        if !bin_member(&v, &u[i]) {
            for j in 0..annx.len() {
                let t = annx[j].2 as usize;
                if !rheaders[t].contains("segment") {
                    let name = rheaders[t].after("|").between("|", "|");
                    if rheaders[t].contains("UTR") && u[i] == name {
                        to_delete[j] = true;
                    }
                }
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Log alignments.

    if verbose {
        fwriteln!(log, "\nALIGNMENTS SIX\n");
        for i in 0..annx.len() {
            print_alignx(log, &annx[i], refdata);
        }
    }

    // In light of the previous calculation, see if one V is aligned much better
    // than another V.  This is done by looking for simple indel events.
    // Probably will have to be generalized.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    let mut vs = Vec::<(usize, usize)>::new();
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if !rheaders[t].contains("V-REGION") {
            continue;
        }
        vs.push((t, i));
    }
    vs.sort_unstable();
    //                     len parts errs  index
    let mut score = Vec::<(i32, usize, usize, usize)>::new();
    let mut j = 0;
    let mut nonsimple = false;
    let mut have_split = false;
    let max_indel = 27;
    let min_len_gain = 100;
    while j < vs.len() {
        let k = next_diff1_2(&vs, j as i32) as usize;
        if k - j == 1 {
            score.push((annx[j].1, k - j, annx[j].4.len(), vs[j].1));
        } else if k - j == 2 {
            let (i1, i2) = (vs[j].1, vs[j + 1].1);
            let (a1, a2) = (&annx[i1], &annx[i2]);
            let mut simple = false;
            let (l1, p1, len1) = (a1.0, a1.3, a1.1);
            let (l2, p2, len2) = (a2.0, a2.3, a2.1);
            if l1 + len1 == l2
                && p1 + len1 < p2
                && (p2 - p1 - len1) % 3 == 0
                && p2 - p1 - len1 <= max_indel
            {
                simple = true;
            }
            if l1 + len1 < l2
                && p1 + len1 == p2
                && (l2 - l1 - len1) % 3 == 0
                && l2 - l1 - len1 <= max_indel
            {
                simple = true;
            }
            if simple {
                have_split = true;
                score.push((len1 + len2, k - j, a1.4.len() + a2.4.len(), vs[j].1));
            } else {
                nonsimple = true;
            }
        } else {
            nonsimple = true;
        }
        j = k;
    }
    if !nonsimple && score.duo() && have_split {
        reverse_sort(&mut score);
        if score[0].0 >= score[1].0 + min_len_gain && score[1].1 == 1 {
            to_delete[score[1].3] = true;
        }
    }
    erase_if(&mut annx, &to_delete);

    // Remove certain subsumed alignments.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        for i2 in 0..annx.len() {
            if i2 == i1 || annx[i1].2 != annx[i2].2 {
                continue;
            }
            let (l1, l2) = (annx[i1].0, annx[i2].0);
            let (len1, len2) = (annx[i1].1, annx[i2].1);
            let (p1, p2) = (annx[i1].3, annx[i2].3);
            if l1 != l2 || p1 != p2 {
                continue;
            }
            if len1 > len2 {
                to_delete[i2] = true;
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // If we see TRBJ1 and not TRBJ2, delete any TRBC2.  And conversely.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    let (mut j1, mut j2) = (false, false);
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if rheaders[t].contains("TRBJ1") {
            j1 = true;
        }
        if rheaders[t].contains("TRBJ2") {
            j2 = true;
        }
    }
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if j1 && !j2 && rheaders[t].contains("TRBC2") {
            to_delete[i] = true;
        }
        if !j1 && j2 && rheaders[t].contains("TRBC1") {
            to_delete[i] = true;
        }
    }
    erase_if(&mut annx, &to_delete);

    // Pick between equally performant Js and likewise for Cs.

    let mut to_delete = vec![false; annx.len()];
    for pass in 0..2 {
        for i1 in 0..annx.len() {
            let t1 = annx[i1].2;
            if pass == 1 {
                if !rheaders[t1 as usize].contains("J-REGION") {
                    continue;
                }
            } else if !rheaders[t1 as usize].contains("C-REGION") {
                continue;
            }
            for i2 in 0..annx.len() {
                let t2 = annx[i2].2;
                if pass == 1 {
                    if !rheaders[t2 as usize].contains("J-REGION") {
                        continue;
                    }
                } else if !rheaders[t2 as usize].contains("C-REGION") {
                    continue;
                }
                let (l1, l2) = (annx[i1].0, annx[i2].0);
                let (len1, len2) = (annx[i1].1, annx[i2].1);
                if l1 != l2 || len1 != len2 {
                    continue;
                }
                let (p1, p2) = (annx[i1].3, annx[i2].3);
                if pass == 1 {
                    if p1 + len1 != refs[t1 as usize].len() as i32 {
                        continue;
                    }
                    if p2 + len2 != refs[t2 as usize].len() as i32 {
                        continue;
                    }
                } else if p1 > 0 || p2 > 0 {
                    continue;
                }
                if annx[i1].4.len() != annx[i2].4.len() {
                    continue;
                }
                if t1 < t2 {
                    to_delete[i2] = true;
                }
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Pick between Cs.

    let mut to_delete = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        let t1 = annx[i1].2;
        if !rheaders[t1 as usize].contains("C-REGION") {
            continue;
        }
        for i2 in 0..annx.len() {
            let t2 = annx[i2].2;
            if !rheaders[t2 as usize].contains("C-REGION") {
                continue;
            }
            let (l1, l2) = (annx[i1].0 as usize, annx[i2].0 as usize);
            let (len1, len2) = (annx[i1].1 as usize, annx[i2].1 as usize);
            // let (p1,p2) = (annx[i1].3,annx[i2].3);
            if l1 + len1 != l2 + len2 {
                continue;
            }
            if l1 + annx[i1].4.len() >= l2 + annx[i2].4.len() {
                continue;
            }
            to_delete[i2] = true;
        }
    }
    erase_if(&mut annx, &to_delete);

    // Again remove UTR annotations that have no matching V annotation.
    // ◼ DANGER with nonstandard references.
    // ◼ Note repetition.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    let (mut u, mut v) = (Vec::<String>::new(), Vec::<String>::new());
    for i in 0..annx.len() {
        let t = annx[i].2 as usize;
        if !rheaders[t].contains("segment") {
            let name = rheaders[t].after("|").between("|", "|");
            if rheaders[t].contains("UTR") {
                u.push(name.to_string());
            }
            if rheaders[t].contains("V-REGION") {
                v.push(name.to_string());
            }
        }
    }
    v.sort();
    for i in 0..u.len() {
        if !bin_member(&v, &u[i]) {
            for j in 0..annx.len() {
                let t = annx[j].2 as usize;
                if !rheaders[t].contains("segment") {
                    let name = rheaders[t].after("|").between("|", "|");
                    if rheaders[t].contains("UTR") && u[i] == name {
                        to_delete[j] = true;
                    }
                }
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Remove some subsumed extended annotations.

    let mut to_delete: Vec<bool> = vec![false; annx.len()];
    for i1 in 0..annx.len() {
        let l1 = annx[i1].0 as usize;
        let len1 = annx[i1].1 as usize;
        for i2 in 0..annx.len() {
            let t2 = annx[i2].2 as usize;
            let l2 = annx[i2].0 as usize;
            let len2 = annx[i2].1 as usize;
            if len2 >= len1 {
                continue;
            }
            if !rheaders[t2].contains("before") && !rheaders[t2].contains("after") {
                continue;
            }
            if l1 <= l2 && l1 + len1 >= l2 + len2 {
                to_delete[i2] = true;
            }
        }
    }
    erase_if(&mut annx, &to_delete);

    // Transform.

    ann.clear();
    for x in annx.iter() {
        ann.push((x.0, x.1, x.2, x.3, x.4.len() as i32));
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// PRINT ANNOTATIONS
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Print annotations, marking any V annotations that are out of frame.

pub fn print_some_annotations(
    refdata: &RefData,
    ann: &[(i32, i32, i32, i32, i32)],
    log: &mut Vec<u8>,
    verbose: bool,
) {
    let refs = &refdata.refs;
    let rheaders = &refdata.rheaders;
    if verbose {
        fwriteln!(log, "");
    }
    let mut vstart = Vec::<i32>::new();
    for l in 0..ann.len() {
        let estart = ann[l].0;
        let t = ann[l].2 as usize;
        let tstart = ann[l].3;
        if tstart == 0 && (rheaders[t].contains("V-REGION") || rheaders[t].contains("L+V")) {
            vstart.push(estart);
        }
    }
    for l in 0..ann.len() {
        let (estart, len) = (ann[l].0, ann[l].1);
        let t = ann[l].2 as usize;
        let tstart = ann[l].3;
        let mis = ann[l].4;
        fwrite!(
            log,
            "{}-{} ==> {}-{} on {} [len={}] (mis={})",
            estart,
            estart + len,
            tstart,
            tstart + len,
            rheaders[t],
            refs[t].len(),
            mis
        );
        if vstart.solo()
            && (rheaders[t].contains("V-REGION") || rheaders[t].contains("L+V"))
            && (estart - vstart[0] - tstart) % 3 != 0
        {
            fwrite!(log, " [SHIFT!]");
        }
        fwriteln!(log, "");
    }
}

pub fn print_annotations(
    b: &DnaString,
    refdata: &RefData,
    log: &mut Vec<u8>,
    allow_improper: bool,
    abut: bool,
    verbose: bool,
) {
    let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
    annotate_seq_core(
        b,
        refdata,
        &mut ann,
        true,
        allow_improper,
        abut,
        log,
        verbose,
    );
    print_some_annotations(refdata, &ann, log, verbose);
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// CDR3
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Given a DNA sequence, return a CDR3 sequence in it (if found), its start
// position on the DNA sequence, and left and right scores (see below).  The CDR3
// sequence is an amino acid sequence having length at least 5, starting with
// a C, and not containing a stop codon.
//
// NOTE THAT THE RIGHT MOTIF OVERLAPS THE CDR3 BY THREE AMINO ACIDS!
//
// In addition, we score the CDR3 and flanking sequences versus left and right
// motifs, and require a minimum score to pass.  These motifs were derived by
// appropriately stacking up V and J segments and looking for high multiplicity
// amino acids at given positions (see jcon.rs).
//
// If more than one CDR3 sequence is found, we first reduce to those having the
// highest score.  Then we choose the ones having the greatest start position.
// Finally we pick the longest motif.
//
// ◼ The interpretation of the tig slice is ugly.  See comments at
// ◼ get_cdr3_using_ann.

pub fn cdr3_motif_left() -> Vec<Vec<u8>> {
    vec![
        b"LQPEDSAVYY".to_vec(),
        b"VEASQTGTYF".to_vec(),
        b"ATSGQASLYL".to_vec(),
    ]
}
pub fn cdr3_motif_right() -> Vec<Vec<u8>> {
    vec![b"LTFG.GTRVTV".to_vec(), b"LIWG.GSKLSI".to_vec()]
}

const CDR3_MIN_LEN: usize = 5;
const LEFT_FLANK_MIN_SCORE: usize = 3;
const RIGHT_FLANK_MIN_SCORE: usize = 4;

pub struct CDR3Annotation {
    pub start_position_on_contig: usize,
    pub aa_seq: Vec<u8>,
    left_flank_score: usize,
    right_flank_score: usize,
}

pub fn get_cdr3(contig: &DnaStringSlice) -> Vec<CDR3Annotation> {
    const MIN_TOTAL_CDR3_SCORE: usize = 10; // about as high as one can go

    let left_motifs = cdr3_motif_left();
    let right_motifs = cdr3_motif_right();

    let contig_seq = contig.to_owned().to_ascii_vec();

    let mut found_cdr3s: Vec<CDR3Annotation> = Vec::new();
    for frame_idx in 0..3 {
        // Go through three frames.

        // Convert the DNA sequence + frame to an amino acid sequence.
        let amino_acid_seq = nucleotide_to_aminoacid_sequence(&contig_seq, frame_idx);
        if amino_acid_seq.len() < 4 {
            continue;
        }

        // Check each position in the AA seq to see if we find a CDR3 there.
        for cdr3_start_pos in
            0..amino_acid_seq.len() - min(amino_acid_seq.len(), (CDR3_MIN_LEN + 3) + 1)
        {
            // The CDR3 has to start with a Cysteine.
            if amino_acid_seq[cdr3_start_pos] == b'C' {
                // Search for the right flank set up start and end positions for the search.
                let first_f = cdr3_start_pos + (CDR3_MIN_LEN - 3);
                let last_f = amino_acid_seq.len() - RIGHT_FLANK_MIN_SCORE;
                for right_motif_start_pos in first_f..last_f {
                    // Don't look further if the remaining part of the amino_acid_seq is to short to find the full right flank motif.
                    if right_motif_start_pos + right_motifs[0].len() > amino_acid_seq.len() {
                        break;
                    }

                    // Match the right flank and calculate the score.
                    let mut right_flank_score = 0;
                    for right_motif_col in 0..right_motifs[0].len() {
                        let mut hit = false;
                        for right_motif_row in &right_motifs {
                            if amino_acid_seq[right_motif_start_pos + right_motif_col]
                                == right_motif_row[right_motif_col]
                            {
                                hit = true;
                            }
                        }
                        if hit {
                            right_flank_score += 1;
                        }
                    }

                    // If right flank score is larger than RIGHT_FLANK_MIN_SCORE, continue and attempt to match left flank,
                    // otherwise continue with next possible CDR3 start position.
                    if right_flank_score >= RIGHT_FLANK_MIN_SCORE {
                        // Check if there is a stop codon in the CDR3.
                        let mut stop_codon = false;
                        for aa in amino_acid_seq
                            .iter()
                            .take(right_motif_start_pos + 2)
                            .skip(cdr3_start_pos + 1)
                        {
                            if aa == &b'*' {
                                stop_codon = true;
                            }
                        }

                        // If there is no stop codon and there is room for the full left motif in the AA seq, match left flank.
                        let ll = left_motifs[0].len();
                        if !stop_codon && cdr3_start_pos >= ll {
                            let mut left_flank_score = 0;
                            for left_motif_col in 0..ll {
                                let mut hit = false;
                                for left_motif_row in &left_motifs {
                                    if amino_acid_seq[cdr3_start_pos - ll + left_motif_col]
                                        == left_motif_row[left_motif_col]
                                    {
                                        hit = true;
                                    }
                                }
                                if hit {
                                    left_flank_score += 1;
                                }
                            }

                            // If the left flank score and total score is above cutoff push the CDR3 to the vec.
                            // ◼ It's possible that the left_flank_score + right_flank_score
                            // ◼ bound should be increased.
                            if left_flank_score >= LEFT_FLANK_MIN_SCORE
                                && left_flank_score + right_flank_score >= MIN_TOTAL_CDR3_SCORE
                            {
                                found_cdr3s.push(CDR3Annotation {
                                    start_position_on_contig: contig.start
                                        + frame_idx
                                        + 3 * cdr3_start_pos,
                                    aa_seq: amino_acid_seq
                                        [cdr3_start_pos..right_motif_start_pos + 2 + 1]
                                        .to_vec(),
                                    left_flank_score: left_flank_score,
                                    right_flank_score: right_flank_score,
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    // Only return cdr3s having the maximum score.
    let max_score = found_cdr3s
        .iter()
        .map(|cdi| cdi.left_flank_score + cdi.right_flank_score)
        .max()
        .unwrap_or(0);
    let to_delete = found_cdr3s
        .iter()
        .map(|cdi| cdi.left_flank_score + cdi.right_flank_score < max_score)
        .collect::<Vec<_>>();
    erase_if(&mut found_cdr3s, &to_delete);
    found_cdr3s.sort_by_key(|cdi| cdi.start_position_on_contig);

    // Prefer later start and prefer longer CDR3.
    if found_cdr3s.len() > 1 {
        // ◼ This is awkward.
        let n = found_cdr3s.len();
        found_cdr3s.swap(0, n - 1);
        found_cdr3s.truncate(1);
    };

    found_cdr3s
}

pub fn print_cdr3(tig: &DnaStringSlice, log: &mut Vec<u8>) {
    let cdr3_anns = get_cdr3(tig);
    for cdr3 in cdr3_anns {
        fwriteln!(
            log,
            "cdr3 = {} at {}, score = {} + {}",
            strme(&cdr3.aa_seq),
            cdr3.start_position_on_contig,
            cdr3.left_flank_score,
            cdr3.right_flank_score
        );
    }
}

// Given annotations of a DNA sequence, return a slice showing where the CDR3
// sequence should live, or a null slice.  This uses some empirically determined
// bounds.
//
// ◼ This seems very unlikely to be optimal.  The value of LOW_RELV_CDR3 was
// ◼ lowered to make BCR work, which suggests that measuring relative to the end
// ◼ of the V segment is not right.

pub fn cdr3_loc<'a>(
    contig: &'a DnaString,
    refdata: &RefData,
    ann: &[(i32, i32, i32, i32, i32)],
) -> DnaStringSlice<'a> {
    // Given the design of this function, the following bound appears to be optimal
    // except possibly for changes less than ten.
    const LOW_RELV_CDR3: isize = -40;
    if ann.is_empty() {
        return contig.slice(0, 0);
    }
    let mut i = ann.len() - 1;
    loop {
        let t = ann[i].2 as usize;
        if !refdata.rheaders[t].contains("segment") && refdata.is_v(t) {
            let (l, p) = (ann[i].0 as isize, ann[i].3 as isize);
            let vstop_on_contig = l + refdata.refs[t].len() as isize - p;
            let mut start = vstop_on_contig + LOW_RELV_CDR3;
            if start < 0 {
                start = 0;
            }
            if start > contig.len() as isize {
                start = contig.len() as isize;
            }
            return contig.slice(start as usize, contig.len());
        }
        if i == 0 {
            return contig.slice(0, 0);
        }
        i -= 1;
    }
}

// Given a DNA sequence and annotations of it, as defined by annotate_seq, find
// CDR3 positions on it, constrained by the annotation.  This uses empirically
// determined bounds relative to that annotation.

pub fn get_cdr3_using_ann(
    tig: &DnaString,
    refdata: &RefData,
    ann: &[(i32, i32, i32, i32, i32)],
) -> Vec<CDR3Annotation> {
    let window = cdr3_loc(tig, refdata, ann);

    // Enlarge the window because get_cdr3 looks for motifs to the left and right
    // of the actual CDR3.
    // ◼ Pretty ugly.  This should really be inside get_cdr3.

    let start = max(0, window.start as isize - cdr3_motif_left()[0].ilen() * 3);
    let mut stop = start + window.length as isize + cdr3_motif_right()[0].ilen() * 3;
    if stop > tig.len() as isize {
        stop = tig.len() as isize;
    }
    if stop < start {
        stop = start;
    }
    let window2 = tig.slice(start as usize, stop as usize);

    // Now find the CDR3.
    get_cdr3(&window2)
}

pub fn print_cdr3_using_ann(
    tig: &DnaString,
    refdata: &RefData,
    ann: &[(i32, i32, i32, i32, i32)],
    log: &mut Vec<u8>,
) {
    let found_cdr3s = get_cdr3_using_ann(tig, refdata, ann);
    for cdr3 in found_cdr3s {
        fwriteln!(
            log,
            "cdr3 = {} at {}, score = {} + {}",
            strme(&cdr3.aa_seq),
            cdr3.start_position_on_contig,
            cdr3.left_flank_score,
            cdr3.right_flank_score
        );
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// ANNOTATION STRUCTURE
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Here we have code that presents annotations.json.

// Coordinates in an AnnotationUnit are zero-based.  The alignment score is computed
// using the following penalties:
// MATCH_SCORE = 2
// MISMATCH_PENALTY = 3
// GAP_OPEN_PENALTY = 5
// EXTEND_PENALTY = 1
// which are copied from cellranger/lib/python/cellranger/vdj/constants.py.

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq)]
pub struct AnnotationFeature {
    pub chain: VdjChain,        // chain type of the reference record, e.g. TRA
    pub display_name: String,   // same as gene_name
    pub feature_id: usize,      // id of reference record
    pub gene_name: String,      // name of reference record e.g. TRAV14-1
    pub region_type: VdjRegion, // region type e.g. L-REGION+V-REGION
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq)]
pub struct AnnotationUnit {
    pub contig_match_start: usize,     // start on contig
    pub contig_match_end: usize,       // stop on contig
    pub annotation_match_start: usize, // start on reference record
    pub annotation_match_end: usize,   // stop on reference record
    pub annotation_length: usize,      // length of reference record
    pub cigar: String,                 // cigar of the alignment
    pub score: i32,                    // score of the alignment
    pub feature: AnnotationFeature,    // feature type
}

impl AnnotationUnit {
    // Given one or two alignment entities as produced by annotate_seq, of the same
    // contig to the same reference sequence, produce an AnnotationUnit.

    pub fn from_annotate_seq(
        b: &DnaString,
        refdata: &RefData,
        ann: &[(i32, i32, i32, i32, i32)],
    ) -> AnnotationUnit {
        // Sanity check the inputs.  Obviously these conditions should be checked
        // before calling, so that they can never fail.

        let na = ann.len();
        assert!(na == 1 || na == 2);
        if ann.len() == 2 {
            assert!(ann[0].2 == ann[1].2);
            assert!(
                (ann[0].0 + ann[0].1 == ann[1].0 && ann[0].3 + ann[0].1 < ann[1].3)
                    || (ann[0].0 + ann[0].1 < ann[1].0 && ann[0].3 + ann[0].1 == ann[1].3)
            );
        }

        // Build a cigar string for a single alignment, having an indel in the case
        // where there are two alignment entities.  This does not show mismatches.

        let mut cig = String::new();
        let left1 = ann[0].0 as usize;
        let len1 = ann[0].1 as usize;
        let right1 = b.len() - left1 - len1;
        if left1 > 0 {
            write!(cig, "{left1}S").unwrap();
        }
        write!(cig, "{len1}M").unwrap();
        if na == 1 && right1 > 0 {
            write!(cig, "{right1}S").unwrap();
        }
        if na == 2 {
            let n1 = ann[1].0 - ann[0].0 - ann[0].1;
            let n2 = ann[1].3 - ann[0].3 - ann[0].1;
            if n1 == 0 {
                write!(cig, "{n2}D").unwrap();
            }
            if n2 == 0 {
                write!(cig, "{n1}I").unwrap();
            }
            let left2 = ann[1].0 as usize;
            let len2 = ann[1].1 as usize;
            let right2 = b.len() - left2 - len2;
            write!(cig, "{len2}M").unwrap();
            if right2 > 0 {
                write!(cig, "{right2}S").unwrap();
            }
        }

        // Test for internal soft clipping, which would be a bug.
        // ◼ This is horrible.  We should have a function validate_cigar_string
        // ◼ that validates a cigar string in its entirety, not just test for one
        // ◼ type of anomaly.

        let mut s_pos = Vec::new();
        let mut char_pos = 0;
        for c in cig.chars() {
            if c.is_ascii_alphabetic() {
                if c == 'S' {
                    s_pos.push(char_pos)
                }
                char_pos += 1;
            }
        }
        for p in &s_pos {
            if *p != 0 && *p != char_pos - 1 {
                panic!("Illegal internal soft clipping in cigar {}", cig);
            }
        }

        // Compute alignment score.

        let mut s = 0_i32;
        let t = ann[0].2 as usize;
        let r = &refdata.refs[t];
        for l in 0..na {
            for i in 0..ann[l].1 {
                if b.get((ann[l].0 + i) as usize) == r.get((ann[l].3 + i) as usize) {
                    s += 2;
                } else {
                    s -= 3;
                }
            }
        }
        if na == 2 {
            let n1 = ann[1].0 - ann[0].0 - ann[0].1;
            let n2 = ann[1].3 - ann[0].3 - ann[0].1;
            let n = max(n1, n2);
            s += 4 + n;
        }

        // Build the rest.

        let types = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
        let mut chain_type = String::new();
        for i in 0..types.len() {
            if refdata.rheaders[t].contains(types[i]) {
                chain_type = types[i].to_string();
                break;
            }
        }
        let v: Vec<&str> = refdata.rheaders[t].split_terminator('|').collect();
        AnnotationUnit {
            contig_match_start: ann[0].0 as usize,
            contig_match_end: (ann[na - 1].0 + ann[na - 1].1) as usize,
            annotation_match_start: ann[0].3 as usize,
            annotation_match_end: (ann[na - 1].3 + ann[na - 1].1) as usize,
            annotation_length: refdata.refs[t].len(),
            cigar: cig,
            score: s,
            feature: AnnotationFeature {
                chain: chain_type.parse().unwrap(),
                display_name: refdata.name[t].clone(),
                feature_id: v[1].force_usize(),
                gene_name: refdata.name[t].clone(),
                region_type: v[3].parse().unwrap(),
            },
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone, Default, PartialEq, Eq)]
pub struct ClonotypeInfo {
    #[serde(default)]
    pub raw_clonotype_id: Option<String>,
    #[serde(default)]
    pub raw_consensus_id: Option<String>,
    #[serde(default)]
    pub exact_subclonotype_id: Option<String>,
}

impl ClonotypeInfo {
    pub fn empty() -> Self {
        ClonotypeInfo::default()
    }
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq)]
pub struct Region {
    pub start: usize,
    pub stop: usize,
    pub nt_seq: String,
    pub aa_seq: String,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq)]
pub struct JunctionSupport {
    pub reads: i32,
    pub umis: i32,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct ContigAnnotation {
    // raw data for the contig
    pub barcode: String,     // the barcode
    pub contig_name: String, // name of the contig
    pub sequence: String,    // nucleotide sequence for contig
    pub quals: String,       // contig quality scores
    pub fraction_of_reads_for_this_barcode_provided_as_input_to_assembly: Option<f64>,

    // contig support
    pub read_count: usize, // number of reads assigned to contig
    pub umi_count: usize,  // number of UMIs assigned to the contig

    // amino acid sequence
    //
    // The start position of the amino acid sequence on the contig is unspecified.
    // ◼ This seems like a flaw.
    pub start_codon_pos: Option<usize>, // start pos on contig of start codon
    pub stop_codon_pos: Option<usize>,  // start pos on contig of stop codon
    pub aa_sequence: Option<String>,    // amino acid sequence
    pub frame: Option<usize>,           // null and never changed (unused field)

    // CDR3
    pub cdr3: Option<String>,      // amino acid sequence for CDR3, or null
    pub cdr3_seq: Option<String>,  // nucleotide sequence for CDR3, or null
    pub cdr3_start: Option<usize>, // start position in bases on contig of CDR3
    pub cdr3_stop: Option<usize>,  // stop position in bases on contig of CDR3

    // CDR* and FWR* regions. There are not filled in here
    // CDR3 follows a different convention here to ensure backwards compatibility
    pub fwr1: Option<Region>,
    pub cdr1: Option<Region>,
    pub fwr2: Option<Region>,
    pub cdr2: Option<Region>,
    pub fwr3: Option<Region>,
    pub fwr4: Option<Region>,

    // annotations
    #[serde(default)]
    pub annotations: Vec<AnnotationUnit>, // the annotations
    pub clonotype: Option<String>, // null, filled in later
    #[serde(default)]
    pub info: ClonotypeInfo, // Empty initially, may be filled in later

    // state of the contig
    pub high_confidence: bool,               // declared high confidence?
    pub validated_umis: Option<Vec<String>>, // validated UMIs
    pub non_validated_umis: Option<Vec<String>>, // non-validated UMIs
    pub invalidated_umis: Option<Vec<String>>, // invalidated UMIs
    pub is_cell: bool,                       // was the barcode declared a cell?
    pub productive: Option<bool>,            // productive?  (null means not full length)
    #[serde(default = "set_true")]
    pub filtered: bool, // true and never changed (unused field)
    /// criteria used to asess productive status
    pub productive_criteria: Option<ContigStatus>,

    pub is_gex_cell: Option<bool>, // Was the barcode declared a cell by Gene expression data, if available
    pub is_asm_cell: Option<bool>, // Was the barcode declared a cell by the VDJ assembler

    pub full_length: Option<bool>, // New field added in CR 4.1. None if the field is not set

    pub junction_support: Option<JunctionSupport>, // New field added in CR 7.2. Coverage of junction region for a good contig
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sample: Option<String>,
}

fn set_true() -> bool {
    true
}

impl ContigAnnotation {
    // Given the alignment entities produced by annotate_seq, produce a
    // ContigAnnotation.  This is done so as to produce at most one V, D, J and C,
    // each.  Pairs of alignment entities that are separated by an indel get
    // collapsed in this process.

    pub fn from_annotate_seq(
        b: &DnaString,                           // the contig
        q: &[u8],                                // qual scores for the contig
        tigname: &str,                           // name of the contig
        refdata: &RefData,                       // reference data
        ann: &[(i32, i32, i32, i32, i32)],       // output of annotate_seq
        nreads: usize,                           // number of reads assigned to contig
        numis: usize,                            // number of umis assigned to contig
        high_confidencex: bool,                  // declared high confidence?
        validated_umis: Option<Vec<String>>,     // validated UMIs
        non_validated_umis: Option<Vec<String>>, // non-validated UMIs
        invalidated_umis: Option<Vec<String>>,   // invalidated UMIs
        is_cellx: bool,                          // was the barcode declared a cell?
        productivex: bool,                       // productive?
        productive_criteria: ContigStatus,       // criteria used to asess productive status
        jsupp: Option<JunctionSupport>,          // num reads, umis supporting junction
    ) -> ContigAnnotation {
        let mut vstart = -1_i32;
        for i in 0..ann.len() {
            let t = ann[i].2 as usize;
            if refdata.is_v(t) && ann[i].3 == 0 {
                vstart = ann[i].0;
            }
        }
        let mut aa = String::new();
        let mut stop = -1_i32;
        let x = b.to_owned().to_ascii_vec();
        if vstart >= 0 {
            let y = nucleotide_to_aminoacid_sequence(&x, vstart as usize);
            aa = stringme(&y);
            for i in 0..y.len() {
                if y[i] == b'*' {
                    stop = vstart + 3 * (i as i32);
                    break;
                }
            }
        }
        let (mut cdr3x, mut cdr3x_dna) = (String::new(), String::new());
        let (mut cdr3x_start, mut cdr3x_stop) = (-1_i32, -1_i32);
        let found_cdr3s = if !refdata.refs.is_empty() {
            get_cdr3_using_ann(b, refdata, ann)
        } else {
            get_cdr3(&b.slice(0, b.len()))
        };
        if !found_cdr3s.is_empty() {
            let cdr3 = found_cdr3s.first().unwrap();
            cdr3x = stringme(&cdr3.aa_seq);
            let start = cdr3.start_position_on_contig;
            for i in start..start + 3 * cdr3x.len() {
                cdr3x_dna.push(x[i] as char);
            }
            cdr3x_start = start as i32;
            cdr3x_stop = (start + 3 * cdr3x.len()) as i32;
        }
        let mut qp = q.to_vec();
        for i in 0..q.len() {
            qp[i] += 33;
        }
        let mut ann = ContigAnnotation {
            barcode: tigname.before("_").to_string(),
            contig_name: tigname.to_string(),
            sequence: b.to_string(),
            quals: stringme(&qp),
            fraction_of_reads_for_this_barcode_provided_as_input_to_assembly: None,
            read_count: nreads,
            umi_count: numis,
            start_codon_pos: match vstart {
                -1 => None,
                _ => Some(vstart as usize),
            },
            stop_codon_pos: match stop {
                -1 => None,
                _ => Some(stop as usize),
            },
            aa_sequence: match vstart {
                -1 => None,
                _ => Some(aa),
            },
            frame: None,
            cdr3: match cdr3x.is_empty() {
                true => None,
                _ => Some(cdr3x.clone()),
            },
            cdr3_seq: match cdr3x.is_empty() {
                true => None,
                _ => Some(cdr3x_dna),
            },
            cdr3_start: match cdr3x.is_empty() {
                true => None,
                _ => Some(cdr3x_start as usize),
            },
            cdr3_stop: match cdr3x.is_empty() {
                true => None,
                _ => Some(cdr3x_stop as usize),
            },
            annotations: make_annotation_units(b, refdata, ann),
            clonotype: None,
            info: ClonotypeInfo::empty(),
            high_confidence: high_confidencex,
            validated_umis,
            non_validated_umis,
            invalidated_umis,
            is_cell: is_cellx,
            productive: Some(productivex),
            productive_criteria: Some(productive_criteria),
            filtered: true,
            junction_support: jsupp,
            // These need to be populated by the assembler explicitly as needed
            is_gex_cell: None,
            is_asm_cell: None,
            full_length: None,

            fwr1: None,
            cdr1: None,
            fwr2: None,
            cdr2: None,
            fwr3: None,
            fwr4: None,

            sample: None,
        };
        ann.full_length = Some(ann.is_full_length());
        ann
    }

    // Produce a ContigAnnotation from a sequence.

    pub fn from_seq(
        b: &DnaString,                           // the contig
        q: &[u8],                                // qual scores for the contig
        tigname: &str,                           // name of the contig
        refdata: &RefData,                       // reference data
        nreads: usize,                           // number of reads assigned to contig
        numis: usize,                            // number of umis assigned to contig
        high_confidence: bool,                   // declared high confidence?
        validated_umis: Option<Vec<String>>,     // validated UMIs
        non_validated_umis: Option<Vec<String>>, // non-validated UMIs
        invalidated_umis: Option<Vec<String>>,   // invalidated UMIs
        is_cell: bool,                           // was the barcode declared a cell?
        jsupp: Option<JunctionSupport>,          // num reads, umis supporting junction
    ) -> ContigAnnotation {
        let mut ann = Vec::<(i32, i32, i32, i32, i32)>::new();
        annotate_seq(b, refdata, &mut ann, true, false, true);
        let (is_productive, productive_criteria) = is_productive_contig(b, refdata, &ann);
        ContigAnnotation::from_annotate_seq(
            b,
            q,
            tigname,
            refdata,
            &ann,
            nreads,
            numis,
            high_confidence,
            validated_umis,
            non_validated_umis,
            invalidated_umis,
            is_cell,
            is_productive,
            productive_criteria,
            jsupp,
        )
    }

    // Output with four space indentation.  Ends with comma and newline.

    pub fn write(&self, out: &mut BufWriter<File>) {
        let buf = Vec::new();
        let formatter = serde_json::ser::PrettyFormatter::with_indent(b"    ");
        let mut ser = serde_json::Serializer::with_formatter(buf, formatter);
        self.serialize(&mut ser).unwrap();
        fwriteln!(out, "{},", String::from_utf8(ser.into_inner()).unwrap());
    }

    // Print.

    pub fn print(&self, log: &mut Vec<u8>) {
        log.append(&mut serde_json::to_vec_pretty(&self).unwrap());
    }

    /// Find annotation unit corresponding to the given region
    pub fn get_region(&self, region: VdjRegion) -> Option<&AnnotationUnit> {
        self.annotations
            .iter()
            .find(|ann_unit| ann_unit.feature.region_type == region)
    }

    /// Find gene name corresponding to the given region if it exists
    pub fn get_gene_name(&self, region: VdjRegion) -> Option<&String> {
        self.get_region(region).map(|unit| &unit.feature.gene_name)
    }

    pub fn is_productive(&self) -> bool {
        self.productive.unwrap_or(false)
    }

    pub fn is_full_length(&self) -> bool {
        self.full_length.unwrap_or_else(|| {
            check_full_length(self.get_region(VdjRegion::V), self.get_region(VdjRegion::J))
        })
    }

    /// The chain corresponding to this contig is defined as the chain type of the V-region
    /// if it exists. This is consistent with the defn used in enclone
    pub fn chain_type(&self) -> Option<VdjChain> {
        self.get_region(VdjRegion::V).map(|unit| unit.feature.chain)
    }
}

fn check_full_length(v_ann: Option<&AnnotationUnit>, j_ann: Option<&AnnotationUnit>) -> bool {
    match (v_ann, j_ann) {
        (Some(v_region), Some(j_region)) => {
            v_region.annotation_match_start == 0
                && j_region.annotation_match_end == j_region.annotation_length
        }
        _ => false,
    }
}

// Given the alignment entities produced by annotate_seq, produce an AnnotationUnit.
// This is done so as to produce at most one V, D, J and C, each.  Pairs of
// alignment entities that are separated by an indel get collapsed in this process.

pub fn make_annotation_units(
    b: &DnaString,
    refdata: &RefData,
    ann: &[(i32, i32, i32, i32, i32)],
) -> Vec<AnnotationUnit> {
    let mut x = Vec::<AnnotationUnit>::new();
    let rtype = &["U", "V", "D", "J", "C"];
    for &rt in rtype {
        let mut locs = Vec::<(usize, usize, usize)>::new();
        let mut j = 0;
        while j < ann.len() {
            let t = ann[j].2 as usize;
            if refdata.segtype[t] != rt {
                j += 1;
                continue;
            }
            let mut entries = 1;
            let mut len = ann[j].1;
            if j < ann.len() - 1
                && ann[j + 1].2 as usize == t
                && ((ann[j].0 + ann[j].1 == ann[j + 1].0 && ann[j].3 + ann[j].1 < ann[j + 1].3)
                    || (ann[j].0 + ann[j].1 < ann[j + 1].0 && ann[j].3 + ann[j].1 == ann[j + 1].3))
            {
                entries = 2;
                len += ann[j + 1].1;
            }
            let mut score = len as usize;
            if refdata.segtype[t] == "V" && ann[j].3 == 0 {
                score += 1_000_000;
            }
            if refdata.segtype[t] == "J" && (ann[j].3 + ann[j].1) as usize == refdata.refs[t].len()
            {
                score += 1_000_000;
            }
            locs.push((score, j, entries));
            j += entries;
        }
        reverse_sort(&mut locs);
        if !locs.is_empty() {
            let (j, entries) = (locs[0].1, locs[0].2);
            let mut annx = Vec::<(i32, i32, i32, i32, i32)>::new();
            for k in j..j + entries {
                annx.push(ann[k]);
            }
            x.push(AnnotationUnit::from_annotate_seq(b, refdata, &annx));
        }
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refx;

    #[test]
    fn test_no_internal_soft_clipping() {
        use refx::RefData;

        let refdata = RefData::from_fasta(String::from(
            "test/inputs/test_no_internal_soft_clipping_ref.fa",
        ));
        // println!("Loaded reference with {} entries", refdata.id.len());

        let contig_seq = DnaString::from_acgt_bytes("AGGAACTGCTCAGTTAGGACCCAGACGGAACCATGGAAGCCCCAGCGCAGCT\
        TCTCTTCCTCCTGCTACTCTGGCTCCCAGATACCACTGGAGAAATAGTGATGACGCAGTCTCCAGCCACCCTGTCTGTGTCTCCAGGGGAAAGAGCC\
        ACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAGCTACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTG\
        CATCCACCAGGGCCACTGGTATCCCAGCCAGGTTCAGTGGCAGTGGGTCTGGGACAGAGTTCACTCTCACCATCAGCAGCCTGCAGTCTGAAGATTT\
        TGCAGTTTATTACTGTCAGCAGTATAATAACTGGCTCATGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAACGAACTGTGGCTGCACCATCT\
        GTCTTCATCTTCCCGCCATCTGATGAGCAGTTGAAATCTGGAACTGCCTCTGTTGTGTGCCTGCTGAATAACTTCTATCCCAGAGAGGCCAAAGTAC\
        AGTGGAAGGTGGATAACGC".as_bytes());

        // Phred quality passed in, convert to raw quality, the ContigAnnotation add the
        // offset back when writing!
        let contig_qual: Vec<u8> = "III]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\
        ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\
        ]]]]]]]]]]]]]]]]IIII!!!IIIIIIIIIIII]]]]]]]]]X]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\
        ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\
        ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\
        ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]\
        ]]".as_bytes().iter().map(|x| x-33).collect();
        let annotation = ContigAnnotation::from_seq(
            &contig_seq,
            &contig_qual,
            "clonotype125_consensus_1",
            &refdata,
            120,
            2,
            true, // high_confidence
            None,
            None,
            None,
            false, // is_cell, should be changed to None
            None,
        );

        // println!("{:#?}", annotation);
        for ann in &annotation.annotations {
            let mut s_pos = Vec::new();
            let mut char_pos = 0;
            for c in ann.cigar.chars() {
                if c.is_ascii_alphabetic() {
                    if c == 'S' {
                        s_pos.push(char_pos)
                    }
                    char_pos += 1;
                }
            }
            if !s_pos.is_empty() {
                println!("Cigar : {:?}", ann.cigar);
                println!("Soft clipping at : {s_pos:?}");
                for p in &s_pos {
                    assert!(*p == 0 || *p == (char_pos - 1))
                }
            }
        }
    }
}
