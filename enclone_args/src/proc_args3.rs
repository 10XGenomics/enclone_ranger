// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// This file contains the two functions proc_xcr and proc_meta.

use enclone_core::defs::{EncloneControl, OriginInfo};
use enclone_core::{expand_integer_ranges, fetch_url, tilde_expand_me};
use io_utils::{dir_list, open_for_read, open_for_write_new, open_userfile_for_read, path_exists};
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::thread;
use std::time;

use string_utils::TextUtils;
use vector_utils::unique_sort;

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

fn expand_analysis_sets(x: &str, ctl: &EncloneControl) -> Result<String, String> {
    let mut tokens = Vec::<String>::new();
    let mut token = String::new();
    for c in x.chars() {
        if c == ',' || c == ':' || c == ';' {
            if !token.is_empty() {
                tokens.push(token.clone());
                token.clear();
            }
            tokens.push(c.to_string());
        } else {
            token.push(c);
        }
    }
    if !token.is_empty() {
        tokens.push(token);
    }
    let mut tokens2 = Vec::<String>::new();
    for token in tokens {
        if let Some(setid) = token.strip_prefix('S') {
            if ctl.gen_opt.internal_run {
                let url = format!("{}/{setid}", ctl.gen_opt.config["sets"]);
                let m = fetch_url(&url)?;
                if m.contains("\"analysis_ids\":[") {
                    let ids = m.between("\"analysis_ids\":[", "]");
                    let mut ids = ids.replace(' ', "");
                    ids = ids.replace('\n', "");
                    let ids = ids.split(',');
                    let mut ids2 = Vec::<&str>::new();

                    // Remove wiped analysis IDs.

                    for id in ids {
                        let url = format!("{}/{id}", ctl.gen_opt.config["ones"]);
                        let m = fetch_url(&url)?;
                        if m.contains("502 Bad Gateway") {
                            return Err(format!(
                                "\nWell, this is sad.  The URL \
                                {url} returned a 502 Bad Gateway \
                                message.  Please try again later or ask someone for help.\n"
                            ));
                        }
                        if !m.contains("\"wiped\"") {
                            ids2.push(id);
                        }
                    }
                    let mut enclone = "~/enclone".to_string();
                    tilde_expand_me(&mut enclone);
                    if path_exists(&enclone) {
                        let mut sets = "~/enclone/sets".to_string();
                        tilde_expand_me(&mut sets);
                        if !path_exists(&sets) {
                            std::fs::create_dir(&sets).unwrap();
                            let mut setid = format!("~/enclone/sets/{setid}");
                            tilde_expand_me(&mut setid);
                            if !path_exists(&setid) {
                                let mut f = open_for_write_new![&setid];
                                let s = format!("{}\n", ids2.iter().format(","));
                                f.write_all(s.as_bytes()).unwrap();
                            }
                        }
                    }

                    // Proceed.

                    for (j, id) in ids2.into_iter().enumerate() {
                        if j > 0 {
                            tokens2.push(",".to_string());
                        }
                        tokens2.push(id.to_string());
                    }
                    continue;
                } else {
                    return Err(format!(
                        "\nIt looks like you've provided an incorrect analysis set ID {setid}.\n"
                    ));
                }
            } else if setid.parse::<usize>().is_ok() {
                let mut set_file = format!("~/enclone/sets/{setid}");
                tilde_expand_me(&mut set_file);
                if path_exists(&set_file) {
                    let mut f = open_for_read![&set_file];
                    let mut s = String::new();
                    f.read_to_string(&mut s).unwrap();
                    s = s.before("\n").to_string();
                    let ids2 = s.split(',');
                    for (j, id) in ids2.enumerate() {
                        if j > 0 {
                            tokens2.push(",".to_string());
                        }
                        tokens2.push(id.to_string());
                    }
                    continue;
                }
            }
        }
        tokens2.push(token);
    }
    let mut y = String::new();
    for t in tokens2 {
        y += t.as_str();
    }
    Ok(y)
}

// Functions to find the path to data.

pub fn get_path_fail(p: &str, ctl: &EncloneControl, source: &str) -> Result<String, String> {
    for x in ctl.gen_opt.pre.iter() {
        let pp = format!("{x}/{p}");
        if path_exists(&pp) {
            return Ok(pp);
        }
    }
    if !path_exists(p) {
        if ctl.gen_opt.pre.is_empty() {
            let path = std::env::current_dir().unwrap();
            return Err(format!(
                "\nIn directory {}, unable to find the path {}.  This came from the {} argument.\n",
                path.display(),
                p,
                source
            ));
        } else {
            let path = std::env::current_dir().unwrap();
            let mut pre_msg =
                "Here are the number of entries in your PRE directories:\n".to_string();
            for x in ctl.gen_opt.pre.iter() {
                let mut count = "(does not exist)".to_string();
                if path_exists(x) {
                    count = dir_list(x).len().to_string();
                }
                writeln!(pre_msg, "{x}: {count}").unwrap();
            }
            return Err(format!(
                "\nIn directory {}, unable to find the\npath {},\n\
                even if prepended by any of the directories \
                in\nPRE={}.\nThis came from the {} argument.\n{}",
                path.display(),
                p,
                ctl.gen_opt.pre.iter().format(","),
                source,
                pre_msg
            ));
        }
    }
    Ok(p.to_string())
}

fn get_path(p: &str, ctl: &EncloneControl, ok: &mut bool) -> String {
    *ok = false;
    for x in ctl.gen_opt.pre.iter() {
        let mut pp = format!("{x}/{p}");
        if pp.starts_with('~') {
            tilde_expand_me(&mut pp);
        }
        if path_exists(&pp) {
            *ok = true;
            return pp;
        }
    }
    let mut pp = p.to_string();
    if pp.starts_with('~') {
        tilde_expand_me(&mut pp);
    }
    *ok = path_exists(&pp);
    pp
}

fn get_path_or_internal_id(
    p: &str,
    ctl: &EncloneControl,
    source: &str,
    spinlock: &Arc<AtomicUsize>,
) -> Result<String, String> {
    if ctl.gen_opt.evil_eye {
        println!("getting path for {p}");
    }
    let mut ok = false;
    let mut pp = get_path(p, ctl, &mut ok);
    if !ok {
        if !ctl.gen_opt.internal_run {
            get_path_fail(&pp, ctl, source)?;
        } else {
            // For internal runs, try much harder.  This is so that internal users can
            // just type an internal numerical id for a dataset and have it always
            // work.  The code that's used here should be placed somewhere else.

            let mut q = p.to_string();
            if q.contains('/') {
                q = q.before("/").to_string();
            }
            if q.parse::<usize>().is_ok() {
                if !ctl.gen_opt.config.contains_key("ones") {
                    let mut msg = "\nSomething is wrong.  This is an internal run, but \
                        the configuration\nvariable \"ones\" is undefined.\n"
                        .to_string();
                    if ctl.gen_opt.config.is_empty() {
                        msg += "In fact, there are no configuration variables.\n";
                    } else {
                        msg += "Here are the configuration variables that are defined:\n\n";
                        for (key, value) in ctl.gen_opt.config.iter() {
                            write!(msg, "{key} = {value}").unwrap();
                        }
                        msg += "\n";
                    }
                    return Err(msg);
                }
                let url = format!("{}/{q}", ctl.gen_opt.config["ones"]);
                // We force single threading around the https access because we observed
                // intermittently very slow access without it.
                while spinlock.load(Ordering::SeqCst) != 0 {}
                spinlock.store(1, Ordering::SeqCst);
                let m = fetch_url(&url)?;
                spinlock.store(0, Ordering::SeqCst);
                if m.contains("502 Bad Gateway") {
                    return Err(format!(
                        "\nWell this is sad.  The URL \
                        {url} yielded a 502 Bad Gateway \
                        message.  Please try again later or ask someone for help.\n"
                    ));
                }
                if m.contains("\"path\":\"") {
                    let path = m.between("\"path\":\"", "\"");
                    if !p.contains('/') {
                        pp = format!("{path}/outs");
                    } else {
                        pp = format!("{path}/{}", p.after("/"));
                    }
                    if !path_exists(&pp) {
                        thread::sleep(time::Duration::from_millis(100));
                        if path_exists(&pp) {
                            return Err(format!(
                                "\nYou are experiencing unstable filesystem access: \
                                100 milliseconds ago, \
                                the path\n\
                                {pp}\nwas not visible, but now it is.  You might consider posting \
                                this problem on an appropriate \
                                the slack channel.\nOr retry again.  enclone is \
                                giving up because \
                                if filesystem access blinks in and out of existence,\n\
                                other more cryptic events are likely to occur.\n"
                            ));
                        } else {
                            return Err(format!(
                                "\nIt looks like you've provided an analysis ID for \
                                which the pipeline outs folder\n{p}\nhas not yet been generated.\n\
                                This path did not exist:\n{pp}\n\n\
                                Here is the stdout:\n{m}\n"
                            ));
                        }
                    }
                } else {
                    return Err(format!(
                        "\nIt looks like you've provided either an incorrect \
                        analysis ID {p} or else one for which\n\
                        the pipeline outs folder has not yet been generated.\n\
                        This URL\n{url}\ndid not provide a path.\n"
                    ));
                }
            } else {
                return Err(format!(
                    "\nAfter searching high and low, your path\n{p}\nfor {source} \
                    cannot be found.\nPlease check its value and also the value \
                    for PRE if you provided that.\n"
                ));
            }
        }
    }
    if !pp.ends_with("/outs") && path_exists(format!("{pp}/outs")) {
        pp = format!("{pp}/outs");
    }
    if ctl.gen_opt.evil_eye {
        println!("path found");
    }
    Ok(pp)
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

// Parse barcode-level information file.

fn parse_bc(mut bc: String, ctl: &mut EncloneControl, call_type: &str) -> Result<(), String> {
    let delimiter;
    let file_type;
    if bc.ends_with(".tsv") {
        delimiter = '\t';
        file_type = "TSV";
    } else {
        delimiter = ',';
        file_type = "CSV";
    }
    let mut origin_for_bc = HashMap::<String, String>::new();
    let mut donor_for_bc = HashMap::<String, String>::new();
    let mut tag = HashMap::<String, String>::new();
    let mut barcode_color = HashMap::<String, String>::new();
    let mut alt_bc_fields = Vec::<(String, HashMap<String, String>)>::new();
    let spinlock: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
    if !bc.is_empty() {
        bc = get_path_or_internal_id(&bc, ctl, call_type, &spinlock)?;
        let f = open_userfile_for_read(&bc);
        let mut first = true;
        let mut fieldnames = Vec::<String>::new();
        let mut barcode_pos = 0;
        let (mut origin_pos, mut donor_pos, mut tag_pos, mut color_pos) = (None, None, None, None);
        let mut to_alt = Vec::<isize>::new();
        for line in f.lines() {
            let s = line.unwrap();
            if first {
                let fields = s.split(delimiter).collect::<Vec<&str>>();
                to_alt = vec![-1_isize; fields.len()];
                if !fields.contains(&"barcode") {
                    let mut origin = "from the bc field used in META";
                    if call_type == "BC" {
                        origin = "from the BC argument";
                    }
                    return Err(format!(
                        "\nThe file\n{bc}\n{origin}\nis missing the barcode field.\n",
                    ));
                }
                for x in fields.iter() {
                    fieldnames.push(x.to_string());
                }
                for i in 0..fields.len() {
                    if fields[i] == "color" {
                        color_pos = Some(i);
                    }
                    if fields[i] == "barcode" {
                        barcode_pos = i;
                    } else if fields[i] == "origin" {
                        origin_pos = Some(i);
                    } else if fields[i] == "donor" {
                        donor_pos = Some(i);
                    } else if fields[i] == "tag" {
                        tag_pos = Some(i);
                    } else {
                        to_alt[i] = alt_bc_fields.len() as isize;
                        alt_bc_fields
                            .push((fields[i].to_string(), HashMap::<String, String>::new()));
                    }
                }
                first = false;
            } else {
                let fields = s.split(delimiter).collect::<Vec<&str>>();
                if fields.len() != fieldnames.len() {
                    let mut origin = "bc in META";
                    if call_type == "BC" {
                        origin = "BC";
                    }
                    return Err(format!(
                        "\nThere is a line\n{}\nin a {} file defined by {}\n\
                         that has {} fields, which isn't right, because the header line \
                         has {} fields.  This is for the file\n{}.\n",
                        s,
                        file_type,
                        origin,
                        fields.len(),
                        fieldnames.len(),
                        bc,
                    ));
                }
                for i in 0..fields.len() {
                    if to_alt[i] >= 0 {
                        alt_bc_fields[to_alt[i] as usize]
                            .1
                            .insert(fields[barcode_pos].to_string(), fields[i].to_string());
                    }
                }
                if !fields[barcode_pos].contains('-') {
                    let mut origin = "bc in META";
                    if call_type == "BC" {
                        origin = "BC";
                    }
                    return Err(format!(
                        "\nThe barcode \"{}\" appears in the file\n{bc}\ndefined \
                         by {origin}.  That doesn't make sense because a barcode\n\
                         should include a hyphen.\n",
                        fields[barcode_pos]
                    ));
                }
                if let Some(origin_pos) = origin_pos {
                    origin_for_bc.insert(
                        fields[barcode_pos].to_string(),
                        fields[origin_pos].to_string(),
                    );
                }
                if let Some(donor_pos) = donor_pos {
                    donor_for_bc.insert(
                        fields[barcode_pos].to_string(),
                        fields[donor_pos].to_string(),
                    );
                }
                if let Some(tag_pos) = tag_pos {
                    tag.insert(fields[barcode_pos].to_string(), fields[tag_pos].to_string());
                }
                if let Some(color_pos) = color_pos {
                    barcode_color.insert(
                        fields[barcode_pos].to_string(),
                        fields[color_pos].to_string(),
                    );
                }
            }
        }
    }
    ctl.origin_info.origin_for_bc.push(origin_for_bc);
    ctl.origin_info.donor_for_bc.push(donor_for_bc);
    ctl.origin_info.tag.push(tag);
    ctl.origin_info.barcode_color.push(barcode_color);
    ctl.origin_info.alt_bc_fields.push(alt_bc_fields);
    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_xcr(
    f: &str,
    gex: &str,
    bc: &str,
    have_gex: bool,
    ctl: &mut EncloneControl,
) -> Result<(), String> {
    ctl.origin_info = OriginInfo::default();
    if ((ctl.gen_opt.tcr || ctl.gen_opt.tcrgd) && f.starts_with("BCR="))
        || ((ctl.gen_opt.bcr || ctl.gen_opt.tcr) && f.starts_with("TCRGD="))
        || ((ctl.gen_opt.bcr || ctl.gen_opt.tcrgd) && f.starts_with("TCR="))
    {
        return Err("\nOnly one of TCR, BCR, or TCRGD can be specified.\n".to_string());
    }
    ctl.gen_opt.tcr = f.starts_with("TCR=");
    ctl.gen_opt.tcrgd = f.starts_with("TCRGD=");
    ctl.gen_opt.bcr = f.starts_with("BCR=");
    let val = if ctl.gen_opt.tcr {
        f.after("TCR=")
    } else if ctl.gen_opt.bcr {
        f.after("BCR=")
    } else if ctl.gen_opt.tcrgd {
        f.after("TCRGD=")
    } else {
        f
    };
    if val.is_empty() {
        return Err(format!(
            "\nYou can't write {f} with no value on the right hand side.\n\
            Perhaps you need to remove some white space from your command line.\n"
        ));
    }
    let val = expand_integer_ranges(val);
    let val = expand_analysis_sets(&val, ctl)?;
    let donor_groups = if ctl.gen_opt.cellranger {
        vec![&val[..]]
    } else {
        val.split(';').collect::<Vec<&str>>()
    };
    let mut gex2 = expand_integer_ranges(gex);
    gex2 = expand_analysis_sets(&gex2, ctl)?;
    let donor_groups_gex = if ctl.gen_opt.cellranger {
        vec![&gex2[..]]
    } else {
        gex2.split(';').collect::<Vec<&str>>()
    };
    let donor_groups_bc = bc.split(';').collect::<Vec<&str>>();
    let xcr = if ctl.gen_opt.bcr {
        "BCR"
    } else if ctl.gen_opt.tcrgd {
        "TCRGD"
    } else {
        "TCR"
    };
    if have_gex && donor_groups_gex.len() != donor_groups.len() {
        return Err(format!(
            "\nThere are {} {} donor groups and {} GEX donor groups, so \
             the {} and GEX arguments do not exactly mirror each \
             other's structure.\n",
            xcr,
            donor_groups.len(),
            donor_groups_gex.len(),
            xcr
        ));
    }
    if !bc.is_empty() && donor_groups_bc.len() != donor_groups.len() {
        return Err(format!(
            "\nThe {xcr} and BC arguments do not exactly mirror each \
             other's structure.\n"
        ));
    }

    for (id, d) in donor_groups.iter().enumerate() {
        let origin_groups = if ctl.gen_opt.cellranger {
            vec![&d[..]]
        } else {
            (*d).split(':').collect::<Vec<&str>>()
        };
        let mut origin_groups_gex = Vec::<&str>::new();
        if have_gex {
            if ctl.gen_opt.cellranger {
                origin_groups_gex = vec![donor_groups_gex[id]];
            } else {
                origin_groups_gex = donor_groups_gex[id].split(':').collect::<Vec<&str>>();
            }
            if origin_groups_gex.len() != origin_groups.len() {
                return Err(format!(
                    "\nFor donor {}, there are {} {} origin groups and {} GEX origin groups, so \
                     the {} and GEX arguments do not exactly mirror each \
                     other's structure.\n",
                    id + 1,
                    xcr,
                    origin_groups.len(),
                    origin_groups_gex.len(),
                    xcr
                ));
            }
        }
        let mut origin_groups_bc = Vec::<&str>::new();
        if !bc.is_empty() {
            origin_groups_bc = donor_groups_bc[id].split(':').collect::<Vec<&str>>();
            if origin_groups_bc.len() != origin_groups.len() {
                return Err(format!(
                    "\nThe {xcr} and BC arguments do not exactly mirror each \
                     other's structure.\n"
                ));
            }
        }
        for (is, s) in origin_groups.iter().enumerate() {
            let mut datasets = if ctl.gen_opt.cellranger {
                vec![&s[..]]
            } else {
                (*s).split(',').collect::<Vec<&str>>()
            };
            for ds in datasets.iter_mut() {
                if ds.ends_with('/') {
                    *ds = ds.rev_before("/");
                }
            }
            let datasets_gex: Vec<&str>;
            let mut datasets_bc = Vec::<&str>::new();
            if have_gex {
                if ctl.gen_opt.cellranger {
                    datasets_gex = vec![origin_groups_gex[is]];
                } else {
                    datasets_gex = origin_groups_gex[is].split(',').collect::<Vec<&str>>();
                }
                if datasets_gex.len() != datasets.len() {
                    return Err(format!(
                        "\nSee {} {} datasets and {} GEX datasets, so \
                         the {} and GEX arguments do not exactly mirror each \
                         other's structure.\n",
                        xcr,
                        datasets.len(),
                        datasets_gex.len(),
                        xcr
                    ));
                }
            }
            if !bc.is_empty() {
                datasets_bc = origin_groups_bc[is].split(',').collect::<Vec<&str>>();
                if datasets_bc.len() != datasets.len() {
                    return Err(format!(
                        "\nThe {xcr} and BC arguments do not exactly mirror each \
                         other's structure.\n"
                    ));
                }
            }
            for (ix, x) in datasets.iter().enumerate() {
                ctl.origin_info.color.push("".to_string());
                ctl.origin_info.tag.push(HashMap::<String, String>::new());
                let donor_name = format!("d{}", id + 1);
                let origin_name = format!("s{}", is + 1);
                ctl.origin_info.donor_id.push(donor_name);
                ctl.origin_info.origin_id.push(origin_name);
                let mut dataset_name = (*x).to_string();
                if dataset_name.contains('/') {
                    dataset_name = dataset_name.rev_after("/").to_string();
                }
                ctl.origin_info.descrips.push(dataset_name.clone());
                ctl.origin_info.dataset_id.push(dataset_name.clone());

                // Now work on the BC path.

                let mut bcx = String::new();
                if !bc.is_empty() {
                    bcx = datasets_bc[ix].to_string();
                }
                parse_bc(bcx, ctl, "BC")?;
            }
        }
    }

    // Get paths.  This will need to change when cellranger switches to multi.  This code is
    // parallelized because this code can indirectly make many calls to path_exists, and the wall
    // clock time for these can add up.  There should be a way to do this that does not involve
    // multithreading.

    let source = if f.contains('=') { f.before("=") } else { f };
    let mut results = Vec::<(String, String, bool, String)>::new();
    for (id, d) in donor_groups.iter().enumerate() {
        let origin_groups = (*d).split(':').collect::<Vec<&str>>();
        let mut origin_groups_gex = Vec::<&str>::new();
        if have_gex {
            origin_groups_gex = donor_groups_gex[id].split(':').collect::<Vec<&str>>();
        }
        for (is, s) in origin_groups.iter().enumerate() {
            let datasets = (*s).split(',').collect::<Vec<&str>>();
            let mut datasets_gex = Vec::<&str>::new();
            if have_gex {
                datasets_gex = origin_groups_gex[is].split(',').collect::<Vec<&str>>();
            }
            for (ix, x) in datasets.iter().enumerate() {
                let p = (*x).to_string();
                let mut pg = String::new();
                if have_gex {
                    pg = datasets_gex[ix].to_string();
                }
                results.push((p, pg, false, String::new()));
            }
        }
    }

    let spinlock: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
    results.par_iter_mut().for_each(|res| {
        let (p, pg) = (&mut res.0, &mut res.1);
        let resx = get_path_or_internal_id(p, ctl, source, &spinlock);
        match resx {
            Err(resx) => res.3 = resx,
            Ok(resx) => {
                *p = resx;
                if ctl.gen_opt.bcr && path_exists(format!("{p}/vdj_b")) {
                    *p = format!("{p}/vdj_b");
                }
                if ctl.gen_opt.bcr && path_exists(format!("{p}/multi/vdj_b")) {
                    *p = format!("{p}/multi/vdj_b");
                }
                if ctl.gen_opt.tcr && path_exists(format!("{p}/vdj_t")) {
                    *p = format!("{p}/vdj_t");
                }
                if ctl.gen_opt.tcr && path_exists(format!("{p}/multi/vdj_t")) {
                    *p = format!("{p}/multi/vdj_t");
                }
                if have_gex {
                    let resx = get_path_or_internal_id(pg, ctl, "GEX", &spinlock);
                    match resx {
                        Err(resx) => res.3 = resx,
                        Ok(resx) => {
                            *pg = resx;
                            if path_exists(format!("{pg}/count")) {
                                *pg = format!("{pg}/count");
                            }
                            if path_exists(format!("{pg}/count_pd")) {
                                *pg = format!("{pg}/count_pd");
                            }
                        }
                    }
                }
            }
        }
    });
    for result in results {
        if !result.3.is_empty() {
            return Err(result.3);
        }
        ctl.origin_info.dataset_path.push(result.0);
        ctl.origin_info.gex_path.push(result.1);
    }

    Ok(())
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn proc_meta_core(lines: &[String], ctl: &mut EncloneControl) -> Result<(), String> {
    let mut fields = Vec::<String>::new();
    let mut donors = Vec::<String>::new();
    for (count, s) in lines.iter().enumerate() {
        if count == 0 {
            fields.extend(s.split(',').map(str::to_string));
            let mut fields_sorted = fields.clone();
            unique_sort(&mut fields_sorted);
            if fields_sorted.len() < fields.len() {
                return Err(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has duplicate field names\nin its first line.\n"
                        .to_string(),
                );
            }
            let allowed_fields = vec![
                "bc".to_string(),
                "bcr".to_string(),
                "donor".to_string(),
                "gex".to_string(),
                "origin".to_string(),
                "tcr".to_string(),
                "tcrgd".to_string(),
                "color".to_string(),
            ];
            for x in fields.iter() {
                if !allowed_fields.contains(x) {
                    return Err(format!(
                        "\nThe CSV file that you specified using the META or METAX argument \
                         has an illegal field name ({x}) in its first line.\n"
                    ));
                }
            }
            ctl.gen_opt.tcr = fields.contains(&"tcr".to_string());
            ctl.gen_opt.tcrgd = fields.contains(&"tcrgd".to_string());
            ctl.gen_opt.bcr = fields.contains(&"bcr".to_string());
            if !ctl.gen_opt.tcr && !ctl.gen_opt.bcr && !ctl.gen_opt.tcrgd {
                return Err(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has neither the field tcr, tcrgd, or bcr in its first line.\n"
                        .to_string(),
                );
            }
            if ctl.gen_opt.tcr && ctl.gen_opt.bcr {
                return Err(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has both the fields tcr and bcr in its first line.\n"
                        .to_string(),
                );
            }
            if ctl.gen_opt.tcr && ctl.gen_opt.tcrgd {
                return Err(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has both the fields tcr and tcrgd in its first line.\n"
                        .to_string(),
                );
            }
            if ctl.gen_opt.bcr && ctl.gen_opt.tcrgd {
                return Err(
                    "\nThe CSV file that you specified using the META or METAX argument \
                     has both the fields tcrgd and bcr in its first line.\n"
                        .to_string(),
                );
            }
        } else if !s.starts_with('#') && !s.is_empty() {
            let val = s.split(',').collect::<Vec<&str>>();
            if val.len() != fields.len() {
                return Err(format!(
                    "\nMETA or METAX file line {} has a different number of fields than the \
                     first line of the file.\n",
                    count + 1
                ));
            }
            let mut path = String::new();
            let mut abbr = String::new();
            let mut gpath = String::new();
            let mut origin = "s1".to_string();
            let mut donor = "d1".to_string();
            let mut color = "".to_string();
            let mut bc = "".to_string();
            for i in 0..fields.len() {
                let x = &fields[i];
                let mut y = val[i].to_string();
                if y.starts_with('"') && y.ends_with('"') {
                    y = y.after("\"").rev_before("\"").to_string();
                }
                if *x == "tcr" || *x == "bcr" || *x == "tcrgd" {
                    if y.contains(':') {
                        path = y.after(":").to_string();
                        abbr = y.before(":").to_string();
                    } else {
                        path = y.to_string();
                        if path.contains('/') {
                            abbr = path.rev_after("/").to_string();
                        } else {
                            abbr = path.clone();
                        }
                    }
                } else if *x == "gex" {
                    gpath = y.to_string();
                } else if *x == "origin" {
                    origin = y.to_string();
                } else if *x == "donor" {
                    donor = y.to_string();
                } else if *x == "color" {
                    color = y.to_string();
                } else if *x == "bc" && !y.is_empty() {
                    bc = y.to_string();
                }
            }

            // Parse bc and finish up.

            parse_bc(bc.clone(), ctl, "META")?;
            let current_ref = false;
            let spinlock: Arc<AtomicUsize> = Arc::new(AtomicUsize::new(0));
            path = get_path_or_internal_id(&path, ctl, "META", &spinlock)?;
            if ctl.gen_opt.bcr && path_exists(format!("{path}/vdj_b")) {
                path = format!("{path}/vdj_b");
            }
            if ctl.gen_opt.bcr && path_exists(format!("{path}/multi/vdj_b")) {
                path = format!("{path}/multi/vdj_b");
            }
            if ctl.gen_opt.tcr && path_exists(format!("{path}/vdj_t")) {
                path = format!("{path}/vdj_t");
            }
            if ctl.gen_opt.tcr && path_exists(format!("{path}/multi/vdj_t")) {
                path = format!("{path}/multi/vdj_t");
            }
            if ctl.gen_opt.tcrgd && path_exists(format!("{path}/vdj_t_gd")) {
                path = format!("{path}/vdj_t_gd");
            }
            if ctl.gen_opt.tcrgd && path_exists(format!("{path}/multi/vdj_t_gd")) {
                path = format!("{path}/multi/vdj_t_gd");
            }
            if !gpath.is_empty() {
                gpath = get_path_or_internal_id(&gpath, ctl, "META", &spinlock)?;
                if path_exists(format!("{gpath}/count")) {
                    gpath = format!("{gpath}/count");
                }
                if path_exists(format!("{gpath}/count_pd")) {
                    gpath = format!("{gpath}/count_pd");
                }
            }
            if current_ref {
                ctl.gen_opt.current_ref = true;
            }
            let dp = donors
                .iter()
                .enumerate()
                .filter_map(|(j, dj)| if donor == *dj { Some(j) } else { None })
                .next();
            if dp.is_none() {
                donors.push(donor.clone());
            }
            ctl.origin_info.descrips.push(abbr.clone());
            ctl.origin_info.dataset_path.push(path);
            ctl.origin_info.gex_path.push(gpath);
            ctl.origin_info.dataset_id.push(abbr);
            ctl.origin_info.donor_id.push(donor);
            ctl.origin_info.origin_id.push(origin);
            ctl.origin_info.color.push(color);
        }
    }
    Ok(())
}

pub fn proc_meta(v: &[String], ctl: &mut EncloneControl) -> Result<(), String> {
    let mut lines_all = Vec::<Vec<String>>::new();
    for f in v.iter() {
        if !path_exists(f) {
            return Err(format!(
                "\nCan't find the file {f} referenced by your META argument.\n"
            ));
        }
        let fx = File::open(f);
        if fx.is_err() {
            return Err(format!(
                "\nProblem with META: unable to read from the file\n\
                 \"{f}\".\nPlease check that that path makes sense and that you have read \
                 permission for it.\n"
            ));
        }
        let f = BufReader::new(fx.unwrap());
        let mut lines = Vec::<String>::new();
        for line in f.lines() {
            let s = line.unwrap();
            lines.push(s);
        }
        lines_all.push(lines);
    }
    let mut lines = Vec::<String>::new();
    for j in 0..lines_all.len() {
        if lines_all[j].is_empty() || lines_all[j][0] != lines_all[0][0] {
            return Err(
                "\nMETA files having different header lines have been specified.\n".to_string(),
            );
        }
        if j == 0 {
            lines.push(lines_all[0][0].clone());
        }
        for k in 1..lines_all[j].len() {
            lines.push(lines_all[j][k].clone());
        }
    }
    proc_meta_core(&lines, ctl)
}
