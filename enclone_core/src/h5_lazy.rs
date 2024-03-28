//! Provide a single type allowing for eager or lazy loading of h5 data.
//! This replaced some sprawling code that was special-cased to u32; for
//! simplicity it is left with explicit types rather than being generic.
use std::ops::Range;

use anyhow::Result;
use hdf5::{Dataset, Reader};

pub struct LazyH5Array<'a> {
    reader: Reader<'a>,
    data: Option<Vec<u32>>,
}

impl<'a> LazyH5Array<'a> {
    /// Initialize a lazy H5 array.
    /// If eager_load is true, the data is immediately loaded into memory rather
    /// than on-demand.
    pub fn new(dataset: &'a Dataset, eager_load: bool) -> Result<Self> {
        let reader = dataset.as_reader();
        let data = if eager_load {
            Some(reader.read_raw()?)
        } else {
            None
        };
        Ok(Self { reader, data })
    }

    /// Return a range of the array as a Vec.
    /// This is inefficient but reflects the implementation of the code this was
    /// extracted from.
    pub fn get_range(&self, range: Range<usize>) -> Result<Vec<u32>> {
        match &self.data {
            Some(d) => Ok(d[range].to_vec()),
            None => Ok(self.reader.read_slice(range)?.to_vec()),
        }
    }
}
