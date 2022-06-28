// Copyright (c) 2021 10X Genomics, Inc. All rights reserved.

use crate::defs::EncloneControl;
use string_utils::{stringme, TextUtils};

#[derive(Clone, PartialEq)]
pub struct LinearCondition {
    pub coeff: Vec<f64>,     // left hand side (lhs) coefficients
    pub var: Vec<String>,    // left hand side variables (parallel to coefficients)
    pub rhs: f64,            // right hand side; sum of lhs must exceed rhs
    pub sense: &'static str, // le, ge, lt, gt
}

impl LinearCondition {
    pub fn n(&self) -> usize {
        self.coeff.len()
    }

    pub fn new(x: &str) -> Result<LinearCondition, String> {
        let y = x.replace(' ', "");
        let lhs: &str;
        let rhs: &str;
        let sense: &str;
        if y.contains(">=") {
            (lhs, rhs) = y.split_once(">=").unwrap();
            sense = "ge";
        } else if y.contains('≥') {
            (lhs, rhs) = y.split_once('≥').unwrap();
            sense = "ge";
        } else if y.contains('⩾') {
            (lhs, rhs) = y.split_once('⩾').unwrap();
            sense = "ge";
        } else if y.contains("<=") {
            (lhs, rhs) = y.split_once("<=").unwrap();
            sense = "le";
        } else if y.contains('≤') {
            (lhs, rhs) = y.split_once('≤').unwrap();
            sense = "le";
        } else if y.contains('⩽') {
            (lhs, rhs) = y.split_once('⩽').unwrap();
            sense = "le";
        } else if y.contains('<') {
            (lhs, rhs) = y.split_once('<').unwrap();
            sense = "lt";
        } else if y.contains('>') {
            (lhs, rhs) = y.split_once('>').unwrap();
            sense = "gt";
        } else {
            return Err(format!(
                "\nImproperly formatted condition, no inequality symbol, \
                 please type \"enclone help display\": {}.\n",
                x
            ));
        }
        let mut rhs = rhs.replace('E', "e");
        if !rhs.contains('.') && !rhs.contains('e') {
            rhs += ".0";
        }
        if rhs.parse::<f64>().is_err() {
            return Err(format!(
                "\nImproperly formatted condition, right-hand side invalid: {}.\n\
                The right-hand side needs to be a constant.  Please type \
                \"enclone help filter\"\n\
                for more information.\n",
                x
            ));
        }
        let rhs = rhs.force_f64();
        let mut parts = Vec::<String>::new();
        let mut last = 0;
        let lhsx = lhs.as_bytes();
        let mut parens = 0_isize;
        for i in 0..lhsx.len() {
            if i > 0 && parens == 0 && (lhsx[i] == b'+' || lhsx[i] == b'-') {
                if lhsx[last] != b'+' {
                    parts.push(stringme(&lhsx[last..i]));
                } else {
                    parts.push(stringme(&lhsx[last + 1..i]));
                }
                last = i;
            }
            if lhsx[i] == b'(' {
                parens += 1;
            } else if lhsx[i] == b')' {
                parens -= 1;
            }
        }
        let mut coeff = Vec::<f64>::new();
        let mut var = Vec::<String>::new();
        if lhsx[last] != b'+' {
            parts.push(stringme(&lhsx[last..]));
        } else {
            parts.push(stringme(&lhsx[last + 1..]));
        }
        for part in parts.iter_mut() {
            *part = part.replace('(', "");
            *part = part.replace(')', "");
            if part.contains('*') {
                let mut coeffi = part.before("*").to_string();
                let vari = part.after("*");
                if !coeffi.contains('.') && !coeffi.contains('e') {
                    coeffi += ".0";
                }
                if coeffi.parse::<f64>().is_err() {
                    return Err(format!(
                        "\nImproperly formatted condition, coefficient {} is invalid: {}.\n\
                        Please type \"enclone help filter\" for more information.\n",
                        coeffi, x
                    ));
                }
                coeff.push(coeffi.force_f64());
                var.push(vari.to_string());
            } else {
                let mut coeffi = 1.0;
                let mut start = 0;
                if part.starts_with('-') {
                    coeffi = -1.0;
                    start = 1;
                }
                coeff.push(coeffi);
                var.push(part[start..].to_string());
            }
        }
        Ok(LinearCondition {
            coeff,
            var,
            rhs,
            sense,
        })
    }

    pub fn satisfied(&self, val: &[f64]) -> bool {
        let lhs: f64 = self
            .coeff
            .iter()
            .zip(val.iter())
            .map(|(&ci, &vi)| ci * vi)
            .sum();
        if self.sense == "lt" {
            lhs < self.rhs
        } else if self.sense == "gt" {
            lhs > self.rhs
        } else if self.sense == "le" {
            lhs <= self.rhs
        } else {
            lhs >= self.rhs
        }
    }

    pub fn require_valid_variables(&self, _ctl: &EncloneControl) -> Result<(), String> {
        for i in 0..self.var.len() {
            if self.var[i].ends_with("_cell") {
                return Err(format!(
                    "\nThe variable {} should not be used in a linear condition.\n\
                    Please type \"enclone help filter\" for more information.\n",
                    self.var[i]
                ));
            }
        }
        Ok(())
    }
}
