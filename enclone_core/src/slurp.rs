// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.
//
// Slurp in needed data from an h5 file.

use hdf5::types::FixedAscii;
use itertools::Itertools;

pub fn slurp_h5(
    h5_path: &str,
    take_matrix: bool,
    barcodes: &mut Vec<String>,
    features: &mut Vec<String>,
    matrix: &mut Vec<Vec<(i32, i32)>>,
) -> Result<(), String> {
    // Read barcodes from the h5 file.

    let h = hdf5::File::open(h5_path).unwrap();
    let barcode_loc = h.dataset("matrix/barcodes").unwrap();

    let barcodes0: Result<Vec<FixedAscii<18>>, hdf5::Error> = barcode_loc.as_reader().read_raw();
    if barcodes0.is_err() {
        return Err(format!(
            "\nencountered error reading HDF5 file\n{h5_path}\nas follows\n{}\n",
            barcodes0.as_ref().err().unwrap()
        ));
    }
    let barcodes0 = barcodes0.unwrap();

    barcodes.extend(barcodes0.into_iter().map_into());

    // Read features from the h5 file.

    let feature_id_loc = h.dataset("matrix/features/id").unwrap();
    let feature_ids: Vec<FixedAscii<256>> = feature_id_loc.as_reader().read_raw().unwrap();
    let feature_name_loc = h.dataset("matrix/features/name").unwrap();
    let feature_names: Vec<FixedAscii<256>> = feature_name_loc.as_reader().read_raw().unwrap();
    let feature_type_loc = h.dataset("matrix/features/feature_type").unwrap();
    let feature_types: Vec<FixedAscii<256>> = feature_type_loc.as_reader().read_raw().unwrap();
    for i in 0..feature_ids.len() {
        features.push(format!(
            "{}\t{}\t{}",
            feature_ids[i], feature_names[i], feature_types[i]
        ));
    }

    // If appropriate, construct the binary matrix file from the h5 file.

    if take_matrix {
        let data_loc = h.dataset("matrix/data").unwrap();
        let data: Vec<u32> = data_loc.as_reader().read_raw().unwrap();
        let ind_loc = h.dataset("matrix/indices").unwrap();
        let ind: Vec<u32> = ind_loc.as_reader().read_raw().unwrap();
        let ind_ptr_loc = h.dataset("matrix/indptr").unwrap();
        let ind_ptr: Vec<u32> = ind_ptr_loc.as_reader().read_raw().unwrap();
        matrix.resize(barcodes.len(), Vec::new());
        for i in 0..matrix.len() {
            for j in ind_ptr[i]..ind_ptr[i + 1] {
                matrix[i].push((ind[j as usize] as i32, data[j as usize] as i32));
            }
        }
    }
    Ok(())
}
