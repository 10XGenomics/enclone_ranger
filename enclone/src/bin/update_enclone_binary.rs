// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Update the internally public enclone binary.  This is also done by start_release, but sometimes
// one wants to update the binary without making a release.

use enclone_core::defs::get_config;
use io_utils::path_exists;

use std::collections::HashMap;
use std::env;

#[cfg(not(target_os = "windows"))]
use std::os::unix::fs::PermissionsExt;

fn main() {
    let mut config = HashMap::<String, String>::new();
    let mut config_file = String::new();
    for (key, value) in env::vars() {
        if key == "ENCLONE_CONFIG" {
            config_file = value.to_string();
        }
    }
    if get_config(&config_file, &mut config) {
        let bin = &config["enclone_linux_bin"];
        if !path_exists(bin) {
            std::fs::create_dir_all(bin).unwrap();
        }
        let current = format!("{bin}/enclone");
        let last = format!("{bin}/enclone_last");
        if path_exists(&last) {
            std::fs::remove_file(&last).unwrap();
        }
        if path_exists(&current) {
            std::fs::rename(&current, &last).unwrap();
        }
        std::fs::copy("target/debug/enclone", &current).unwrap();
        #[cfg(not(target_os = "windows"))]
        {
            let perms = std::fs::Permissions::from_mode(0o775);
            std::fs::set_permissions(&current, perms).unwrap();
        }
    }
}
