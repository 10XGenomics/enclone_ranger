// Copyright (c) 2020 10X Genomics, Inc. All rights reserved.

// Rather than put this in a build.rs script in the enclone_proto crate to
// auto generate types.rs from types.proto, this is used to update types.rs
// offline.  A unit test ensures that they're in sync.
// This allows dependent crates to avoid having `prost_build` in their
// transitive dependency tree, and also makes fresh builds (e.g. in CI)
// for dependent crates quite a lot faster.

use prost_build::Config;
use std::path::{Path, PathBuf};
use std::process::{exit, Command};

fn main() {
    let manifest_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let out_dir = manifest_dir.join("../enclone_proto/src");
    make_output(manifest_dir.as_path(), out_dir.as_path());
}

fn make_output(manifest_dir: &Path, out_dir: &Path) {
    std::env::set_var("OUT_DIR", out_dir.as_os_str());
    let mut config = Config::new();
    config.type_attribute(".", "#[derive(::serde::Serialize, ::serde::Deserialize)]");
    config
        .compile_protos(&[manifest_dir.join("types.proto")], &[manifest_dir])
        .unwrap();

    let status = Command::new("rustfmt")
        .arg(out_dir.join("enclone.types.rs"))
        .status()
        .expect("failed to execute rustfmt");
    if !status.success() {
        println!("rustfmt did not complete successfully!");
        exit(status.code().unwrap_or_default());
    }
}

#[cfg(test)]
mod test {
    use super::make_output;
    use std::{fs::read_to_string, path::PathBuf};

    // Ensure that the checked-in file is up to date.
    #[test]
    fn check_output_unchanged() {
        let manifest_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
        let out_dir = tempfile::tempdir().unwrap();
        make_output(manifest_dir.as_path(), out_dir.path());
        let current_source =
            read_to_string(manifest_dir.join("enclone.types.rs").as_path()).unwrap();
        let new_source = read_to_string(out_dir.path().join("enclone.types.rs").as_path()).unwrap();
        assert_eq!(current_source, new_source);
    }
}
