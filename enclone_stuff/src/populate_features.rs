// Copyright (c) 2021 10x Genomics, Inc. All rights reserved.

// Populate features.

use amino::aa_seq;
use enclone_core::defs::EncloneControl;
use io_utils::fwriteln;
use std::fmt::Write as _;
use std::io::Write;
use string_utils::{stringme, strme};
use tables::print_tabular_vbox;
use vdj_ann::refx::RefData;
use vdj_ann::vdj_features::{cdr1_start, cdr2_start, fr1_start, fr2_start, fr3_start};

pub fn populate_features(
    ctl: &EncloneControl,
    refdata: &RefData,
    broken: &[bool],
    fr1_starts: &mut Vec<usize>,
    fr2_starts: &mut Vec<Option<usize>>,
    fr3_starts: &mut Vec<Option<usize>>,
    cdr1_starts: &mut Vec<Option<usize>>,
    cdr2_starts: &mut Vec<Option<usize>>,
    log: &mut Vec<u8>,
) -> Result<(), String> {
    *fr1_starts = vec![0; refdata.refs.len()];
    *fr2_starts = vec![None; refdata.refs.len()];
    *fr3_starts = vec![None; refdata.refs.len()];
    *cdr1_starts = vec![None; refdata.refs.len()];
    *cdr2_starts = vec![None; refdata.refs.len()];
    let mut msg = String::new();
    for i in 0..refdata.refs.len() {
        if refdata.is_v(i) {
            if broken[i] && ctl.gen_opt.require_unbroken_ok {
                continue;
            }
            let aa = aa_seq(&refdata.refs[i].to_ascii_vec(), 0);
            let rtype = refdata.rtype[i];
            let chain_type = if rtype == 0 {
                "IGH"
            } else if rtype == 1 {
                "IGK"
            } else if rtype == 2 {
                "IGL"
            } else if rtype == 3 {
                "TRA"
            } else if rtype == 4 {
                "TRB"
            } else {
                continue;
            };
            let fs1 = fr1_start(&aa, chain_type);
            fr1_starts[i] = 3 * fs1;
            let fs2 = fr2_start(&aa, chain_type, false);
            if let Some(fs2) = fs2 {
                fr2_starts[i] = Some(3 * fs2);
            } else if ctl.gen_opt.require_unbroken_ok {
                msg += "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR2 start \
                    could not be computed\nfor this reference sequence:";
                let seq = refdata.refs[i].to_ascii_vec();
                write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                    .expect("formatting to string");
            }
            let fs3 = fr3_start(&aa, chain_type, false);
            if let Some(fs3) = fs3 {
                fr3_starts[i] = Some(3 * fs3);
            } else if ctl.gen_opt.require_unbroken_ok {
                msg += "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR3 start \
                    could not be computed\nfor this reference sequence:";
                let seq = refdata.refs[i].to_ascii_vec();
                write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                    .expect("formatting to string");
            }
            let cs1 = cdr1_start(&aa, chain_type, false);
            if let Some(cs1) = cs1 {
                cdr1_starts[i] = Some(3 * cs1);
                if let Some(fs2) = fs2 {
                    if cs1 > fs2 && ctl.gen_opt.require_unbroken_ok {
                        msg +=
                            "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR1 start \
                        exceeds the FWR2 start for this reference sequence:";
                        let seq = refdata.refs[i].to_ascii_vec();
                        write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                            .expect("formatting to string");
                    }
                }
            } else if ctl.gen_opt.require_unbroken_ok {
                msg += "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR1 start \
                    could not be computed\nfor this reference sequence:\n";
                let seq = refdata.refs[i].to_ascii_vec();
                write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                    .expect("formatting to string");
            }
            let cs2 = cdr2_start(&aa, chain_type, false);
            if let Some(cs2) = cs2 {
                cdr2_starts[i] = Some(3 * cs2);
                if let Some(fs3) = fs3 {
                    if ctl.gen_opt.require_unbroken_ok && cs2 > fs3 {
                        msg +=
                            "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR2 start \
                        exceeds the FWR3 start for this reference sequence:";
                        let seq = refdata.refs[i].to_ascii_vec();
                        write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                            .expect("formatting to string");
                    }
                }
            } else if ctl.gen_opt.require_unbroken_ok {
                msg += "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the CDR2 start \
                    could not be computed\nfor this reference sequence:";
                let seq = refdata.refs[i].to_ascii_vec();
                write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                    .expect("formatting to string");
            }
            if let Some(cs1) = cs1 {
                if fs1 > cs1 && ctl.gen_opt.require_unbroken_ok {
                    msg += "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR1 start \
                    exceeds the CDR1 start for this reference sequence:\n";
                    let seq = refdata.refs[i].to_ascii_vec();
                    write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                        .expect("formatting to string");
                }
            }
            if let (Some(cs2), Some(fs2)) = (cs2, fs2) {
                if fs2 > cs2 && ctl.gen_opt.require_unbroken_ok {
                    msg += "\nYou supplied the argument REQUIRE_UNBROKEN_OK, but the FWR2 start \
                    exceeds the CDR2 start for this reference sequence:";
                    let seq = refdata.refs[i].to_ascii_vec();
                    write!(msg, ">{}\n{}\n", refdata.rheaders_orig[i], strme(&seq))
                        .expect("formatting to string");
                }
            }
        }
    }
    if !msg.is_empty() {
        return Err(msg);
    }

    // Report on broken reference sequences.  This comes after the json loading because possibly
    // the user supplied the wrong reference, so there is no value in criticizing the reference
    // in that case.

    if !log.is_empty() && !ctl.gen_opt.cellranger && !ctl.gen_opt.accept_broken {
        let mut log = Vec::<u8>::new();
        fwriteln!(
            log,
            "\nSome errors were detected in the reference sequences supplied to enclone.\n\
            Please see comments at end for what you can do about this.\n\n",
        );
        fwriteln!(log,
"🌼  Dear user, some defects were detected in the reference sequences supplied to enclone.   🌼\n\
 🌼  Some of these defects may be small.  Generally they are associated with V segments that 🌼\n\
 🌼  are frameshifted or truncated, or with C segments that have an extra base at the        🌼\n\
 🌼  beginning.  We are letting you know about this because they could result in             🌼\n\
 🌼  misannotation.                                                                          🌼\n"
        );

        let rows = [
            vec![
                "You can make enclone ignore these defects by adding the additional argument"
                    .to_string(),
            ],
            vec![
                "ACCEPT_BROKEN to the enclone command line.  Or you can obtain the same"
                    .to_string(),
            ],
            vec![
                "behavior by defining the environment variable ENCLONE_ACCEPT_BROKEN.".to_string(),
            ],
        ];
        let mut log = stringme(&log);
        print_tabular_vbox(&mut log, &rows, 2, b"l".as_ref(), false, true);
        let mut log = log.as_bytes().to_vec();
        fwriteln!(log, "");
        fwriteln!(
            log,
            "This is probably OK, but if your sample is human or mouse, you may wish to either:\n\
        • rerun cellranger using the cleaned up reference sequences that come prepackaged with \
          it\n  (noting that your might have used an older, less clean version of that)\n\
        • or add the argument BUILT_IN to enclone, which will force use of the built-in \
        reference\n  \
        sequences.  This will be a bit slower because all the contigs will need to be\n  \
        reannotated.  If you're using mouse, you'll also need to add the argument MOUSE.\n"
        );
        return Err(stringme(&log));
    }
    Ok(())
}
