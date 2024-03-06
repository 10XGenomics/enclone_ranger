// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use self::annotate::{annotate_seq, get_cdr3_using_ann, print_some_annotations, Annotation};
use self::refx::RefData;
use self::transcript::is_productive_contig;
use debruijn::dna_string::DnaString;
use enclone_core::barcode_fate::BarcodeFate;
use enclone_core::defs::{EncloneControl, OriginInfo, TigData};
use io_utils::{open_maybe_compressed, path_exists};
use martian_filetypes::json_file::{Json, LazyJsonReader};
use martian_filetypes::LazyRead;
use rand::Rng;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write;
use std::io::BufReader;
use string_utils::{stringme, strme, TextUtils};
use vdj_ann::annotate::ContigAnnotation;
use vdj_ann::{annotate, refx, transcript};
use vdj_types::{VdjChain, VdjRegion};
use vector_utils::{bin_position, erase_if, unique_sort};

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn json_error(json: Option<&str>, internal_run: bool, msg: &str) -> String {
    let mut msgx =
        "There is something wrong with the contig annotations in the cellranger output file"
            .to_string();
    if let Some(json) = json {
        write!(msgx, "\n{json}.").unwrap();
    } else {
        msgx += ".";
    }
    if internal_run {
        writeln!(msgx, "\n\npossibly relevant internal data: {msg}").unwrap();

        msgx += "\n\nATTENTION INTERNAL 10X USERS!\n\
            Quite possibly you are using data from a cellranger run carried out using a \
            version\n\
            between 3.1 and 4.0.  For certain of these versions, it is necessary to add the\n\
            argument CURRENT_REF to your command line.  If that doesn't work, \
            please see below.\n";
    }
    msgx += "\n\nHere is what you should do:\n\n\
        1. If you used cellranger version ≥ 4.0, the problem is very likely\n\
        that the directory outs/vdj_reference was not retained, so enclone\n\
        didn't see it, and had to guess what the reference sequence was.\n\
        Fix this and everything should be fine.\n\n\
        2. If you used cellranger version 3.1, then you need to add a command-line\n\
        argument REF=<vdj_reference_fasta_file_name>, or if you already did that,\n\
        make sure it is the *same* as that which you gave cellranger.\n\n\
        3. If you used cellranger version < 3.1 (the only other possibility), then\n\
        you have options:\n\
        • rerun cellranger using the current version\n\
        • or provide an argument REF= as above and RE to force reannotation\n\
        • or provide the argument BUILT_IN to use the current reference and force\n  \
            reannotation (and MOUSE if you used mouse); only works with human and mouse.\n\n\
        Note that one way to get the error is to specify TCR when you meant BCR, or the\n\
        other way.\n\n\
        If you're stuck, please write to us at enclone@10xgenomics.com.\n";

    msgx
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

#[derive(Default)]
struct JsonParseResult {
    vdj_cell: Option<String>,
    gex_cell: Option<String>,
    gex_cells_specified: bool,
    tig: Option<TigData>,
}

fn process_json_annotation(
    ann: ContigAnnotation,
    json: &str,
    accept_inconsistent: bool,
    origin_info: &OriginInfo,
    dataset_index: usize,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    reannotate: bool,
    ctl: &EncloneControl,
) -> Result<JsonParseResult, String> {
    let mut res: JsonParseResult = Default::default();

    // Get cell status.  Sometime after CR 4.0 was released, and before 4.1 was released,
    // we added new fields is_asm_cell and is_gex_cell to the json file.  The value of
    // is_asm_cell is the original determination of "cell" in the VDJ pipeline, whereas the
    // value of is_gex_cell is that for the GEX pipeline.
    let mut is_cell = ann.is_cell;
    if ann.is_asm_cell.is_some_and(|is_asm_cell| is_asm_cell) {
        is_cell = true;
    }

    if let Some(is_gex_cell) = ann.is_gex_cell {
        res.gex_cells_specified = true;
        if is_gex_cell {
            res.gex_cell = Some(ann.barcode.clone());
        }
    }

    if !ctl.gen_opt.ncell && !is_cell {
        return Ok(res);
    }
    if is_cell {
        res.vdj_cell = Some(ann.barcode.clone());
    }

    // Proceed.

    if !ctl.gen_opt.reprod && !ann.productive.unwrap_or(false) {
        return Ok(res);
    }
    if !ctl.gen_opt.reprod && !ctl.gen_opt.ncell && !ann.high_confidence {
        return Ok(res);
    }
    let mut left = false;
    let (mut v_ref_id, mut j_ref_id) = (1000000, 0);
    let mut d_ref_id: Option<usize> = None;
    let mut c_ref_id = None;
    let mut chain_type = String::new();
    let mut u_ref_id = None;
    let (mut tig_start, mut tig_stop) = (-1_isize, -1_isize);
    let mut v_stop = 0;
    let mut v_stop_ref = 0;
    let mut d_start = None;
    let mut j_start = 0;
    let mut j_start_ref = 0;
    let mut c_start = None;
    let mut annv = Vec::<Annotation>::new();
    let mut cdr3_aa: String;
    let mut cdr3_dna: String;
    let mut cdr3_start: usize;

    let frac_reads_used = ann
        .fraction_of_reads_for_this_barcode_provided_as_input_to_assembly
        .map(|f| (f * 1_000_000.0).round() as u32);

    // Reannotate.
    if reannotate || ctl.gen_opt.reprod {
        let x = DnaString::from_dna_string(&ann.sequence);
        let mut ann1 = Vec::<Annotation>::new();
        annotate_seq(&x, refdata, &mut ann1, true, false, true);

        // If there are multiple V segment alignments, possibly reduce to just one.

        let mut ann2 = Vec::<Annotation>::new();
        let mut j = 0;
        while j < ann1.len() {
            // let t = ann1[j].ref_tig as usize;
            let mut k = j + 1;
            while k < ann1.len() {
                if refdata.segtype[ann1[k].ref_tig as usize]
                    != refdata.segtype[ann1[j].ref_tig as usize]
                {
                    break;
                }
                k += 1;
            }
            if refdata.segtype[ann1[j].ref_tig as usize] == "V" && k - j > 1 {
                let mut entries = 1;
                if j < ann1.len() - 1
                    && ann1[j + 1].ref_tig as usize == ann1[j].ref_tig as usize
                    && ((ann1[j].seq_start + ann1[j].match_len == ann1[j + 1].seq_start
                        && ann1[j].ref_start + ann1[j].match_len < ann1[j + 1].ref_start)
                        || (ann1[j].seq_start + ann1[j].match_len < ann1[j + 1].seq_start
                            && ann1[j].ref_start + ann1[j].match_len == ann1[j + 1].ref_start))
                {
                    entries = 2;
                }
                ann2.extend(&ann1[j..j + entries]);
            } else {
                ann2.extend(&ann1[j..k]);
            }
            j = k;
        }
        ann1 = ann2;

        // Proceed.

        if ctl.gen_opt.trace_barcode == ann.barcode {
            let mut log = Vec::<u8>::new();
            print_some_annotations(refdata, &ann1, &mut log, false);
            print!("\n{}", strme(&log));
        }
        if ctl.gen_opt.trace_barcode == ann.barcode {
            if !is_productive_contig(&x, refdata, &ann1).0 {
                println!("invalid");
                return Ok(res);
            }
        } else if !is_productive_contig(&x, refdata, &ann1).0 {
            return Ok(res);
        }
        let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
        get_cdr3_using_ann(&x, refdata, &ann1, &mut cdr3);
        cdr3_aa = stringme(&cdr3[0].1);
        cdr3_start = cdr3[0].0;
        cdr3_dna = x
            .slice(cdr3_start, cdr3_start + 3 * cdr3_aa.len())
            .to_string();
        let mut seen_j = false;
        for anni in ann1 {
            // let t = anni.ref_tig as usize;
            if refdata.is_u(anni.ref_tig as usize) {
                u_ref_id = Some(anni.ref_tig as usize);
            } else if refdata.is_v(anni.ref_tig as usize) && !seen_j {
                v_ref_id = anni.ref_tig as usize;
                annv.push(anni);
                chain_type = refdata.name[anni.ref_tig as usize][0..3].to_string();
                if chain_type == *"IGH"
                    || chain_type == *"TRB"
                    || (chain_type == *"TRD" && ctl.gen_opt.gamma_delta)
                {
                    left = true;
                }
                if anni.ref_start == 0 {
                    tig_start = anni.seq_start as isize;
                    if tig_start > cdr3_start as isize {
                        panic!(
                            "Something is wrong with the CDR3 start for this contig:\n\n{}.",
                            ann.sequence
                        );
                    }
                    cdr3_start -= tig_start as usize;
                }
                v_stop = (anni.seq_start + anni.match_len) as usize;
                v_stop_ref = (anni.ref_start + anni.match_len) as usize;
            } else if refdata.is_d(anni.ref_tig as usize) {
                d_start = Some(anni.seq_start as usize);
                d_ref_id = Some(anni.ref_tig as usize);
            } else if refdata.is_j(anni.ref_tig as usize) {
                j_ref_id = anni.ref_tig as usize;
                tig_stop = (anni.seq_start + anni.match_len) as isize;
                j_start = anni.seq_start as usize;
                j_start_ref = anni.ref_start as usize;
                seen_j = true;
            } else if refdata.is_c(anni.ref_tig as usize) {
                c_ref_id = Some(anni.ref_tig as usize);
                c_start = Some(anni.seq_start as usize);
            }
        }
        for i in (0..annv.len()).rev() {
            annv[i].seq_start -= annv[0].seq_start;
        }
    } else {
        // Use annotations from json file.

        cdr3_aa = ann.cdr3.unwrap();
        cdr3_dna = ann.cdr3_seq.unwrap();
        cdr3_start = ann.cdr3_start.unwrap();
        let annotations = ann.annotations;
        if annotations.is_empty() {
            return Err(format!(
                "\nThe file\n{json}\ndoes not contain annotations.  To use enclone with it, \
                    please specify the argument BUILT_IN\nto force use of the internal \
                    reference and recompute annotations.\n"
            ));
        }
        let mut cigarv = String::new(); // cigar for V segment
        for a in annotations {
            let region_type = a.feature.region_type;
            let feature_id = a.feature.feature_id;
            if !to_ref_index.contains_key(&feature_id) {
                continue;
            }
            let feature_idx = to_ref_index[&feature_id];
            let ref_start = a.annotation_match_start;
            if region_type == VdjRegion::V {
                v_stop = a.contig_match_end;
                v_stop_ref = a.annotation_match_end;
            }
            let gene_name = a.feature.gene_name;
            if refdata.name[feature_idx] != gene_name && !accept_inconsistent {
                return Err(format!(
                    "\nThere is an inconsistency between the reference \
                     file used to create the Cell Ranger output files in\n{}\nand the \
                     reference that enclone is using.\n\nFor example, the feature \
                     numbered {} is\nthe gene {} in one and the gene {} in the other.\n\n\
                     As far as we know, this type of error can only occur with Cell Ranger \
                     versions before 4.0.\n\n\
                     If this is mouse data, please use the argument MOUSE, and that may \
                     solve the problem.\n\n\
                     If this is human or mouse data, and you are OK with using the current \
                     built-in reference that\nenclone has, \
                     you can instead add the argument BUILT_IN to the command line.  This \
                     forces\nrecomputation of annotations and may be somewhat slower.\n\n\
                     A solution that should always work is to supply\n\
                     REF=vdj_reference_fasta_filename as an argument to enclone.\n",
                    json.rev_before("/"),
                    feature_id,
                    gene_name,
                    refdata.name[feature_idx]
                ));
            }
            if region_type == VdjRegion::V && ref_start == 0 {
                let chain = a.feature.chain;
                chain_type = chain.to_string();
                tig_start = a.contig_match_start as isize;
                cdr3_start -= tig_start as usize;
                if chain == VdjChain::IGH
                    || chain == VdjChain::TRB
                    || (chain == VdjChain::TRD && ctl.gen_opt.gamma_delta)
                {
                    left = true;
                }
                v_ref_id = feature_idx;
                cigarv = a.cigar;
            } else {
                // also check for IG chain?????????????????????????????????????????
                let ref_stop = a.annotation_match_end;
                let ref_len = a.annotation_length;
                if region_type == VdjRegion::J && ref_stop == ref_len {
                    tig_stop = a.contig_match_end as isize;
                    j_ref_id = feature_idx;
                    j_start = a.contig_match_start;
                    j_start_ref = a.annotation_match_start;
                }
                if region_type == VdjRegion::UTR {
                    u_ref_id = Some(feature_idx);
                }
                if region_type == VdjRegion::D {
                    d_start = Some(a.contig_match_start);
                    d_ref_id = Some(feature_idx);
                }
                if region_type == VdjRegion::C {
                    c_ref_id = Some(feature_idx);
                    c_start = Some(a.contig_match_start);
                }
            }
        }
        if v_ref_id == 1000000 {
            return Ok(res);
        }

        // Compute annv from cigarv.  We don't compute the mismatch entry.

        let mut cg = Vec::<Vec<u8>>::new(); // pieces of cigar string
        let mut piece = Vec::<u8>::new();
        for c in cigarv.chars() {
            piece.push(c as u8);
            if c.is_ascii_alphabetic() {
                cg.push(piece.clone());
                piece.clear();
            }
        }
        let t = v_ref_id as i32;
        let (mut len1, mut len2) = (0, 0);
        let (mut ins, mut del) = (0, 0);
        for cgi in cg {
            let x = strme(&cgi[0..cgi.len() - 1]).force_i32();
            if cgi[cgi.len() - 1] == b'M' {
                if len1 == 0 {
                    len1 = x;
                } else if len2 == 0 {
                    len2 = x;
                } else {
                    // probably can't happen
                    len1 = 0;
                    len2 = 0;
                    break;
                }
            }
            if cgi[cgi.len() - 1] == b'I' {
                ins = x;
            }
            if cgi[cgi.len() - 1] == b'D' {
                del = x;
            }
        }
        annv.push(Annotation {
            seq_start: 0,
            match_len: len1,
            ref_tig: t,
            ref_start: 0,
            mismatches: 0,
        });
        if ins > 0 && ins % 3 == 0 && del == 0 && len2 > 0 {
            let start = len1 + ins;
            annv.push(Annotation {
                seq_start: start,
                match_len: len2,
                ref_tig: t,
                ref_start: len1,
                mismatches: 0,
            })
        } else if del > 0 && del % 3 == 0 && ins == 0 && len2 > 0 {
            annv.push(Annotation {
                seq_start: len1,
                match_len: len2,
                ref_tig: t,
                ref_start: len1 + del,
                mismatches: 0,
            })
        }
        let rt = &refdata.refs[v_ref_id];
        if annv.len() == 2 && annv[0].match_len as usize > rt.len() {
            let msg = format!("annv[0].1 = {}, rt.len() = {}", annv[0].match_len, rt.len());
            return Err(json_error(None, ctl.gen_opt.internal_run, &msg));
        }

        // Check to see if the CDR3 sequence has changed.  This could happen if the cellranger
        // version for all_contig_annotations.json used an older version of the CDR3 calculation
        // than is used in the current version of enclone.  This could result in internal
        // inconsistencies, leading to an assert somewhere downstream.

        let mut cdr3 = Vec::<(usize, Vec<u8>, usize, usize)>::new();
        let x = DnaString::from_dna_string(&ann.sequence);
        get_cdr3_using_ann(&x, refdata, &annv, &mut cdr3);
        if cdr3.is_empty() {
            return Ok(res);
        }
        let cdr3_aa_alt = stringme(&cdr3[0].1);
        if cdr3_aa != cdr3_aa_alt {
            // This is particularly pathological and rare:

            if tig_start as usize > cdr3[0].0 {
                return Ok(res);
            }

            // Define start.

            cdr3_start = cdr3[0].0 - tig_start as usize;

            // Define cdr3.

            cdr3_aa = cdr3_aa_alt;
            cdr3_dna = x
                .slice(cdr3_start, cdr3_start + 3 * cdr3_aa.len())
                .to_string();
        }
    }

    // Test for two very rare conditions where the CDR3 is busted.  This could be confusing to
    // users if they hit one of these.
    // Case 1: seen on 47680, barcode CGCCAAGTCCATGAAC-1.
    // Case 2: seen on 99640, barcode CAGTAACCATGTCGAT-1.
    // It is not known if these correspond to bugs in cellranger that were subsequently fixed.

    if cdr3_aa.contains('*') {
        return Ok(res);
    }
    if cdr3_start + 3 * cdr3_aa.len() > tig_stop as usize - tig_start as usize {
        return Ok(res);
    }

    // Keep going.

    if tig_start < 0 || tig_stop < 0 {
        let msg = format!("tig_start = {tig_start}, tig_stop = {tig_stop}");
        return Err(json_error(Some(json), ctl.gen_opt.internal_run, &msg));
    }
    let (tig_start, tig_stop) = (tig_start as usize, tig_stop as usize);
    let mut quals = ann.quals.as_bytes().to_vec();
    assert_eq!(ann.sequence.len(), ann.quals.as_bytes().len());
    let seq = &ann.sequence[tig_start..tig_stop].to_string();
    for qual in &mut quals {
        *qual -= 33_u8;
    }
    let full_quals = quals;
    let quals = full_quals[tig_start..tig_stop].to_vec();
    let umi_count = ann.umi_count;
    let read_count = ann.read_count;
    let origin = origin_info.origin_for_bc[dataset_index]
        .get(&ann.barcode)
        .or_else(|| {
            // the way we use s1 here is flaky
            if !origin_info.origin_id[dataset_index].is_empty()
                && (origin_info.origin_id[dataset_index] != *"s1"
                    || origin_info.origin_for_bc[dataset_index].is_empty())
            {
                Some(&origin_info.origin_id[dataset_index])
            } else {
                None
            }
        });
    let donor = origin_info.donor_for_bc[dataset_index]
        .get(&ann.barcode)
        .or_else(|| {
            // the way we use d1 here is flaky
            if !origin_info.origin_id[dataset_index].is_empty()
                && (origin_info.donor_id[dataset_index] != *"d1"
                    || origin_info.donor_for_bc[dataset_index].is_empty())
            {
                Some(&origin_info.donor_id[dataset_index])
            } else {
                None
            }
        });
    let tag = origin_info.tag[dataset_index].get(&ann.barcode);
    let mut origin_index = None;
    let mut donor_index = None;
    let mut tag_index = None;
    if let Some(origin) = origin {
        origin_index = Some(bin_position(&origin_info.origin_list, origin) as usize);
        if let Some(donor) = donor {
            donor_index = Some(bin_position(&origin_info.donor_list, donor) as usize);
        }
    }
    if let Some(tag) = tag {
        tag_index = Some(bin_position(&origin_info.tag_list, tag) as usize);
    }

    res.tig = Some(TigData {
        cdr3_dna,
        len: seq.len(),
        v_start: tig_start,
        v_stop,
        v_stop_ref,
        d_start,
        j_start,
        j_start_ref,
        j_stop: tig_stop,
        c_start,
        full_seq: ann.sequence.as_bytes().to_vec(),
        v_ref_id,
        d_ref_id,
        j_ref_id,
        c_ref_id,
        u_ref_id,
        fr1_start: 0,
        cdr1_start: None,
        fr2_start: None,
        cdr2_start: None,
        fr3_start: None,
        cdr3_aa,
        cdr3_start,
        quals,
        full_quals,
        barcode: ann.barcode,
        tigname: ann.contig_name,
        left,
        dataset_index,
        origin_index,
        donor_index,
        tag_index,
        umi_count,
        read_count,
        chain_type,
        annv,
        validated_umis: ann.validated_umis,
        non_validated_umis: ann.non_validated_umis,
        invalidated_umis: ann.invalidated_umis,
        frac_reads_used,
    });
    Ok(res)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse the JSON annotations file.
//
// Tracking contigs using bc_cdr3_aa; could improve later.
//
// This section requires 3.1.  If you want to avoid that, do something to make tig_start
// and tig_stop always nonnegative.  Or use the RE option.

#[derive(Default)]
struct ReadJsonResult {
    vdj_cells: Vec<String>,
    gex_cells: Vec<String>,
    gex_cells_specified: bool,
    tig_bc: Vec<Vec<TigData>>,
}

fn read_json(
    accept_inconsistent: bool,
    origin_info: &OriginInfo,
    dataset_index: usize,
    json: &String,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
    reannotate: bool,
    ctl: &EncloneControl,
) -> Result<ReadJsonResult, String> {
    let mut jsonx = json.clone();
    if !path_exists(json) {
        jsonx = format!("{json}.lz4");
    }
    if jsonx.contains('/') {
        let p = jsonx.rev_before("/");
        if !path_exists(p) {
            return Err(format!(
                "\nThere should be a directory\n\
                 \"{p}\"\n\
                 but it does not exist.  Please check how you have specified the\n\
                 input files to enclone, including the PRE argument.\n"
            ));
        }
    }
    if !path_exists(&jsonx) {
        return Err(format!(
            "\nThe path\n\
             \"{jsonx}\"\n\
             does not exist.  Please check how you have specified the\n\
             input files to enclone, including the PRE argument.\n"
        ));
    }

    let mut tigs = Vec::new();
    let mut vdj_cells = Vec::new();
    let mut gex_cells = Vec::new();
    let mut gex_cells_specified = false;

    let reader: LazyJsonReader<ContigAnnotation, Json, _> =
        LazyJsonReader::with_reader(BufReader::new(open_maybe_compressed(&jsonx)))
            .map_err(|err| format!("{err:#?}"))?;

    for entry in reader.into_iter() {
        let result = process_json_annotation(
            entry.map_err(|err| err.to_string())?,
            json,
            accept_inconsistent,
            origin_info,
            dataset_index,
            refdata,
            to_ref_index,
            reannotate,
            ctl,
        )?;
        if let Some(tig) = result.tig {
            tigs.push(tig);
        }
        if let Some(c) = result.vdj_cell {
            vdj_cells.push(c);
        }
        if let Some(c) = result.gex_cell {
            gex_cells.push(c);
        }
        if result.gex_cells_specified {
            gex_cells_specified = true;
        }
    }
    unique_sort(&mut gex_cells);
    let mut tig_bc = Vec::<Vec<TigData>>::new();
    let mut r = 0;
    while r < tigs.len() {
        let mut s = r + 1;
        while s < tigs.len() {
            if tigs[s].barcode != tigs[r].barcode {
                break;
            }
            s += 1;
        }

        // For now we require at most four contigs (but we don't yet merge foursies).

        if s - r <= 4 || ctl.clono_filt_opt_def.nmax {
            let mut bc_tigs = tigs[r..s].to_vec();
            bc_tigs.sort();
            tig_bc.push(bc_tigs);
        }
        r = s;
    }
    unique_sort(&mut vdj_cells);

    // Subsample.

    if ctl.gen_opt.subsample >= 0.0 {
        let mut rng = rand::thread_rng();
        let mut to_delete1 = vec![false; tig_bc.len()];
        let mut to_delete2 = vec![false; vdj_cells.len()];
        let mut to_delete3 = vec![false; gex_cells.len()];
        for (bc, del) in tig_bc.iter().zip(to_delete1.iter_mut()) {
            let y: f64 = rng.gen();
            if y < 1.0 - ctl.gen_opt.subsample {
                *del = true;
                let bc = &bc[0].barcode;
                let p = bin_position(&vdj_cells, bc);
                if p >= 0 {
                    to_delete2[p as usize] = true;
                }
                let p = bin_position(&gex_cells, bc);
                if p >= 0 {
                    to_delete3[p as usize] = true;
                }
            }
        }
        erase_if(&mut tig_bc, &to_delete1);
        erase_if(&mut vdj_cells, &to_delete2);
        erase_if(&mut gex_cells, &to_delete3);
    }

    // Done.

    Ok(ReadJsonResult {
        vdj_cells,
        gex_cells,
        gex_cells_specified,
        tig_bc,
    })
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub struct Annotations {
    pub vdj_cells: Vec<Vec<String>>,
    pub gex_cells: Vec<Vec<String>>,
    pub gex_cells_specified: Vec<bool>,
    pub tig_bc: Vec<Vec<TigData>>,
    pub fate: Vec<HashMap<String, BarcodeFate>>,
}

pub fn parse_json_annotations_files(
    ctl: &EncloneControl,
    refdata: &RefData,
    to_ref_index: &HashMap<usize, usize>,
) -> Result<Annotations, String> {
    // Note: only tracking truncated seq and quals initially
    let ann = if !ctl.gen_opt.cellranger {
        "all_contig_annotations.json"
    } else {
        "contig_annotations.json"
    };
    let results = ctl
        .origin_info
        .dataset_path
        .par_iter()
        .enumerate()
        .map(|(li, dataset_path)| {
            let json = format!("{dataset_path}/{ann}");
            let json_lz4 = format!("{dataset_path}/{ann}.lz4");
            if !path_exists(&json) && !path_exists(&json_lz4) {
                return Err(format!("\ncan't find {json} or {json_lz4}\n"));
            }
            read_json(
                ctl.gen_opt.accept_inconsistent,
                &ctl.origin_info,
                li,
                &json,
                refdata,
                to_ref_index,
                ctl.gen_opt.reannotate,
                ctl,
            )
            .map(|r| (li, r))
        })
        .collect::<Result<Vec<_>, String>>()?;

    let mut ann = Annotations {
        tig_bc: Default::default(),
        vdj_cells: Default::default(),
        gex_cells: Default::default(),
        gex_cells_specified: Default::default(),
        fate: vec![HashMap::<String, BarcodeFate>::new(); ctl.origin_info.n()],
    };

    for (i, result) in results {
        let cells = &result.vdj_cells;
        let mut found = vec![false; cells.len()];
        let tigs = &result.tig_bc;
        for tig in tigs {
            let p = bin_position(cells, &tig[0].barcode);
            if p >= 0 {
                found[p as usize] = true;
            }
        }
        for j in 0..found.len() {
            if !found[j] {
                ann.fate[i].insert(cells[j].clone(), BarcodeFate::NonProductive);
            }
        }

        ann.tig_bc.extend(result.tig_bc.into_iter());
        ann.vdj_cells.push(result.vdj_cells);
        ann.gex_cells.push(result.gex_cells);
        ann.gex_cells_specified.push(result.gex_cells_specified);
    }
    Ok(ann)
}
