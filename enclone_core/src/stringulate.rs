// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

// Representation of Vec<String> objects as Strings.  We do this under that assumption that the
// strings in the vector do not contain double escapes, which seems like a reasonable assumption
// in typical circumstances.  Then the representation of v as a string is
// (double escape) v[0] (double escape) v[1] (double escape) ... (double escape).
//
// For this to be reversible, one has to know the length of v.  In fact, the way we use this
// representation is as follows.  Let have a data structure, say Widget, and we convert it to
// a Vec<String>, in such a way that the vector starts with
// "Widget", number of entries in the string, ... .
//
// Second, in practice what we want to represent are vectors of heterogeneous objects, from a
// list of types having functions to_string and from_string, using the above design.  In that
// case we do this:
// 1. Convert each of the objects to strings, except for strings, which we leave intact.
// 2. Concatenate these.

use std::fmt::Display;

use itertools::Itertools;
use string_utils::TextUtils;

const DOUBLE: &str = "";

pub fn flatten_vec_string(v: &[&str]) -> String {
    format!("{DOUBLE}{}{DOUBLE}", v.iter().format(DOUBLE))
}

pub fn unflatten_string(s: &str) -> Vec<&str> {
    let mut chars = s.char_indices();
    // Skip first two characters and the last character.
    let start = chars.nth(2).unwrap().0;
    let end = chars.nth_back(1).unwrap().0;
    let mid = &s[start..end];
    mid.split(&DOUBLE).collect()
}

pub struct HetString {
    pub name: String,
    pub content: String,
}

pub fn unpack_to_het_string(s: &str) -> Vec<HetString> {
    let mut v = Vec::<HetString>::new();
    let fields: Vec<&str> = s.split(DOUBLE).collect();
    let mut i = 0;
    while i < fields.len() {
        if !fields[i].is_empty() {
            v.push(HetString {
                name: "String".to_string(),
                content: fields[i].to_string(),
            });
        }
        if i + 2 < fields.len() {
            let n = fields[i + 2].force_usize();
            v.push(HetString {
                name: fields[i + 1].to_string(),
                content: flatten_vec_string(&fields[i + 1..i + 1 + n]),
            });
            i += n + 1;
        } else {
            break;
        }
    }
    v
}

// Specific implementations, should split off if significantly more are added.

pub struct DescriptionTable {
    pub display_text: String,
    pub spreadsheet_text: String,
}

impl Display for DescriptionTable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}DescriptionTable{}4{}{}{}{}{}",
            DOUBLE, DOUBLE, DOUBLE, self.display_text, DOUBLE, self.spreadsheet_text, DOUBLE,
        )
    }
}

impl DescriptionTable {
    pub fn from_string(x: &str) -> Self {
        let v = unflatten_string(x);
        DescriptionTable {
            display_text: v[2].to_string(),
            spreadsheet_text: v[3].to_string(),
        }
    }
}

#[derive(Default, Clone)]
pub struct FeatureBarcodeAlluvialTable {
    pub id: String,
    pub display_text: String,
    pub spreadsheet_text: String,
}

#[derive(Default, Clone)]
pub struct FeatureBarcodeAlluvialTableSet {
    pub s: Vec<FeatureBarcodeAlluvialTable>,
}

impl Display for FeatureBarcodeAlluvialTableSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}FeatureBarcodeAlluvialTableSet{}{}{}",
            DOUBLE,
            DOUBLE,
            3 * self.s.len() + 2,
            DOUBLE
        )?;
        for s in &self.s {
            write!(
                f,
                "{}{DOUBLE}{}{DOUBLE}{}{DOUBLE}",
                s.id, s.display_text, s.spreadsheet_text
            )?;
        }
        Ok(())
    }
}

impl FeatureBarcodeAlluvialTableSet {
    pub fn from_string(x: &str) -> Self {
        let v = unflatten_string(x);
        let n = v[1].force_usize() / 3;
        let mut s = Vec::with_capacity(n);
        for i in 0..n {
            s.push(FeatureBarcodeAlluvialTable {
                id: v[2 + 3 * i].to_string(),
                display_text: v[2 + 3 * i + 1].to_string(),
                spreadsheet_text: v[2 + 3 * i + 2].to_string(),
            });
        }
        FeatureBarcodeAlluvialTableSet { s }
    }
}

#[derive(Default, Clone)]
pub struct FeatureBarcodeAlluvialReadsTable {
    pub id: String,
    pub display_text: String,
    pub spreadsheet_text: String,
}

#[derive(Default, Clone)]
pub struct FeatureBarcodeAlluvialReadsTableSet {
    pub s: Vec<FeatureBarcodeAlluvialReadsTable>,
}

impl Display for FeatureBarcodeAlluvialReadsTableSet {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}FeatureBarcodeAlluvialReadsTableSet{}{}{}",
            DOUBLE,
            DOUBLE,
            3 * self.s.len() + 2,
            DOUBLE
        )?;
        for s in &self.s {
            write!(
                f,
                "{}{DOUBLE}{}{DOUBLE}{}{DOUBLE}",
                s.id, s.display_text, s.spreadsheet_text
            )?;
        }
        Ok(())
    }
}

impl FeatureBarcodeAlluvialReadsTableSet {
    pub fn from_string(x: &str) -> Self {
        let v = unflatten_string(x);
        let n = v[1].force_usize() / 3;
        let mut s = Vec::with_capacity(n);
        for i in 0..n {
            s.push(FeatureBarcodeAlluvialReadsTable {
                id: v[2 + 3 * i].to_string(),
                display_text: v[2 + 3 * i + 1].to_string(),
                spreadsheet_text: v[2 + 3 * i + 2].to_string(),
            });
        }
        FeatureBarcodeAlluvialReadsTableSet { s }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_unflatten_string() {
        // All ascii.
        let simple_case = "0123456789";
        let simple = unflatten_string(simple_case);
        assert_eq!(simple.len(), 1);
        assert_eq!(simple[0], &simple_case[2..simple_case.len() - 2]);
        // multi-byte characters.
        let complex_case = format!("ƒƒƒ2{DOUBLE}345{DOUBLE}67ƒƒƒ");
        let result = unflatten_string(complex_case.as_str());
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], "ƒ2");
        assert_eq!(result[1], "345");
        assert_eq!(result[2], "67ƒ");
    }

    #[test]
    fn test_flatten_vec_string() {
        let expected = format!("{DOUBLE}ƒ123ƒ{DOUBLE}ƒƒ{DOUBLE}");
        let split = unflatten_string(expected.as_str());
        assert_eq!(expected, flatten_vec_string(&split));
    }
}
