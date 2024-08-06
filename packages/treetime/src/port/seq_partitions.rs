use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::port::seq_sparse::VarPos;
use crate::seq::range::RangeCollection;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct SeqPartition<'a, 'g> {
  gtr: &'g GTR,
  length: usize,
  alphabet: &'a Alphabet,
}

impl<'a, 'g> SeqPartition<'a, 'g> {
  pub fn new(gtr: &'g GTR, length: usize, alphabet: &'a Alphabet) -> Self {
    SeqPartition { gtr, length, alphabet }
  }

  pub fn profile(&self, c: char) -> &Array1<f64> {
    self.alphabet.get_profile(c)
  }

  pub fn code(&self, p: &Array1<f64>) -> char {
    self.alphabet.get_code(p)
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeqInfoLh {
  unknown: RangeCollection,
  gaps: RangeCollection,
  non_char: RangeCollection, // any position that does not evolve according to the substitution model, i.e. gap or N
  variable: BTreeMap<usize, VarPos>, // probability vector for each variable position collecting information from children
  log_lh: f64,
  fixed: BTreeMap<String, Array1<f64>>, // probability vector for the state of fixed positions based on information from children
  fixed_composition: BTreeMap<String, usize>, // this is important for likelihood calculations
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SeqInfoParsimony {
  unknown: RangeCollection,
  gaps: RangeCollection,
  non_char: RangeCollection, // any position that does not evolve according to the substitution model, i.e. gap or N
  variable: BTreeMap<usize, VarPos>, // 0/1 vector for each variable position collecting information from children
  fixed_composition: BTreeMap<String, usize>, // this is important for likelihood calculations
}
