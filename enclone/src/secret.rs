// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Get counts for secreted and membrane proteins.
//
// The basic problem with this is the relevant exon junctions are too far from the 5' end,
// so there is not enough signal to be useful.
//
// Also what is done below appears not to handle the IGHG case correctly.
//
// For IGHG, the boundaries should be
// (A)CH2-(B)CH3-CHS  [secreted]
// (A)CH2-(B)Mx [membrane].

use enclone_core::defs::EncloneControl;
use std::process::Command;
use std::{collections::HashMap, path::Path};
use string_utils::{strme, TextUtils};
use vector_utils::next_diff1_3;

// copied from tenkit2/pack_dna.rs:

pub fn reverse_complement(x: &mut [u8]) {
    x.reverse();
    for xj in x {
        *xj = match *xj {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => *xj,
        }
    }
}

pub fn fetch_secmem(ctl: &mut EncloneControl) -> Result<(), String> {
    // Define the CH3 exon boundaries, and the sequences that could follow it, both in
    // GRCh38 or GRCm38 coordinates.

    let species = &ctl.gen_opt.species;
    let ch3;
    let fol;
    if species == "human" {
        ch3 = [
            ('-', "chr14:105600482-105600805"),
            ('-', "chr14:105840368-105840691"),
            ('-', "chr14:105854918-105855235"),
        ];
        fol = [
            ("TACCTG", "M1"),
            ("GTGAAA", "M2"),
            ("GTGAAG", "M2"),
            ("TGTGAA", "M2?"),
            ("GGCTCT", "M"),
            ("GGAAAC", "S"),
            ("TATGTA", "S"),
            ("GGCCCG", "S"),
            ("GCCCGC", "S"),
            ("GGACAG", "S"),
            ("GGGGTG", "S"),
        ]
        .as_ref();
    } else {
        ch3 = [
            ('-', "chr12:113414273-113414593"),
            ('-', "chr12:113271711-113272031"),
            ('-', "chr12:113421370-113421686"),
        ];
        fol = [
            ("GAGCTAGAC", "M1"),
            ("GAGCTGGAA", "M1"),
            ("GAGGGGGAG", "M1"),
            ("GGCATAGTC", "M1"),
            ("GGGCTAGAC", "M1"),
            ("GGGCTGCAA", "M1"),
            ("GTGAAA", "M2"),
            ("GTGAAG", "M2"),
            ("GAACGTCAA", "M"),
            ("GGAAGAGCC", "S"),
            ("GGCAGACCG", "S"),
            ("GGGCCAGTA", "S"),
            ("GGGCTAGTC", "S"),
            ("GGGTCAGTA", "S"),
            ("GTGAACACC", "S"),
            ("TGAACACCT", "S"),
            ("GAGGTGCAC", "S"),
            ("GCCAGCGCT", "S"),
            ("GGCCAGCGC", "S"),
        ]
        .as_ref();
    }

    // Traverse the datasets.

    for gex_path in ctl.origin_info.gex_path.iter().take(ctl.origin_info.n()) {
        let mut data = Vec::<(String, String, String)>::new(); // (barcode, umi, class)
        let bam = Path::new(gex_path).join("possorted_genome_bam.bam");

        // Traverse the boundaries.

        for ch3i in ch3 {
            // Call samtools.

            let o = Command::new("samtools")
                .arg("view")
                .arg(&bam)
                .arg(ch3i.1)
                .output()
                .expect("failed to execute samtools");

            // Parse the output.

            let o = String::from_utf8(o.stdout).unwrap();
            for line in o.lines() {
                let fields = line.split('\t').collect::<Vec<&str>>();
                let pos = fields[3].force_usize();
                let cigar = fields[5];
                let seq = fields[9];
                let (mut barcode, mut umi) = ("", "");
                for &fj in &fields[11..] {
                    if fj.starts_with("CB:Z:") {
                        barcode = fj.after("CB:Z:");
                    } else if fj.starts_with("UB:Z:") {
                        umi = fj.after("UB:Z:");
                    }
                }
                if barcode.is_empty() {
                    continue;
                }

                // Parse cigar string.

                let mut cg = Vec::<Vec<u8>>::new(); // pieces of cigar string
                let mut piece = Vec::<u8>::new();
                for c in cigar.chars() {
                    piece.push(c as u8);
                    if c.is_ascii_alphabetic() {
                        cg.push(piece.clone());
                        piece.clear();
                    }
                }
                if !piece.is_empty() {
                    cg.push(piece);
                }

                // Determine if the sequence is reaching off the end of the reference interval.
                // This is for the left end in the rc case, and right end otherwise.  The latter
                // case has not been tested.

                let mut ref_pos = pos;
                let mut read_pos = 1;
                let low = ch3i.1.after(":").before("-").force_usize();
                let high = ch3i.1.after("-").force_usize();
                let mut ext = 0;
                let mut ext_seq = Vec::<u8>::new();
                for cgj in cg {
                    let x = cgj[cgj.len() - 1];
                    let n = strme(&cgj[..cgj.len() - 1]).force_usize();
                    if x == b'M' {
                        if ch3i.0 == '-' {
                            if read_pos > 1
                                && ref_pos < high
                                && ref_pos + n > low
                                && read_pos + low > ref_pos + 1
                            {
                                ext = read_pos + low - ref_pos - 1;
                                ext_seq = seq.as_bytes()[0..ext].to_vec();
                                reverse_complement(&mut ext_seq);
                                break;
                            }
                        } else if ref_pos <= high && ref_pos + n > high {
                            ext = ref_pos + n - high;
                            ext_seq = seq.as_bytes()[seq.len() - ext..].to_vec();
                            break;
                        }
                        ref_pos += n;
                        read_pos += n;
                    } else if x == b'N' || x == b'S' || x == b'I' || x == b'D' {
                        ref_pos += n;
                    } else {
                        return Err("\nUnexpected character in cigar string.\n".to_string());
                    }
                }

                // Check if extension long enough.

                if (species == "human" && ext < 6) || (species == "mouse" && ext < 9) {
                    continue;
                }

                // Print.

                let mut class = if species == "human" {
                    strme(&ext_seq[0..6])
                } else {
                    strme(&ext_seq[0..9])
                };
                for &fj in fol {
                    if strme(&ext_seq).starts_with(fj.0) {
                        class = fj.1;
                    }
                }
                data.push((barcode.to_string(), umi.to_string(), class.to_string()));
            }
        }

        // Fill in the map.

        let mut h = HashMap::<String, (usize, usize)>::new();
        data.sort();
        let mut i = 0;
        while i < data.len() {
            let j = next_diff1_3(&data, i as i32) as usize;
            let (mut sec, mut mem) = (0, 0);
            let mut k = i;
            while k < j {
                // let l = next_diff12_3(&data, k as i32) as usize; // crashed loader
                let mut l = k;
                while l < j {
                    if data[l].1 != data[k].1 {
                        break;
                    }
                    l += 1;
                }
                let (mut s, mut m) = (0, 0);
                for dz in &data[k..l] {
                    if dz.2.starts_with('M') {
                        m += 1;
                    } else if dz.2.starts_with('S') {
                        s += 1;
                    }
                }
                if s > 0 && m == 0 {
                    sec += 1;
                } else if s == 0 && m > 0 {
                    mem += 1;
                }
                k = l;
            }
            h.insert(data[i].0.to_string(), (sec, mem));
            i = j;
        }
        ctl.origin_info.secmem.push(h);
    }
    Ok(())
}
