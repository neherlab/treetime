use crate::gtr::gtr::GTR;
use crate::port::seq_sparse::VarPos;
use crate::seq::range::RangeCollection;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct SeqPartition<'g> {
  pub gtr: &'g GTR,
  pub length: usize,
}

impl<'g> SeqPartition<'g> {
  pub fn new(gtr: &'g GTR, length: usize) -> Self {
    SeqPartition { gtr, length }
  }

  pub fn profile(&self, c: char) -> &Array1<f64> {
    self.gtr.alphabet().get_profile(c)
  }

  pub fn code(&self, p: &Array1<f64>) -> char {
    self.gtr.alphabet().get_code(p)
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeqInfoLh {
  pub unknown: RangeCollection,
  pub gaps: RangeCollection,
  pub non_char: RangeCollection, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub variable: BTreeMap<usize, VarPos>, // probability vector for each variable position collecting information from children
  pub log_lh: f64,
  pub fixed: BTreeMap<String, Array1<f64>>, // probability vector for the state of fixed positions based on information from children
  pub fixed_composition: BTreeMap<String, usize>, // this is important for likelihood calculations
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeqInfoParsimony {
  pub unknown: RangeCollection,
  pub gaps: RangeCollection,
  pub non_char: RangeCollection, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub variable: BTreeMap<usize, VarPos>, // 0/1 vector for each variable position collecting information from children
  pub fixed_composition: BTreeMap<String, usize>, // this is important for likelihood calculations
}
