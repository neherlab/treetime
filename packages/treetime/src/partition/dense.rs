use crate::alphabet::alphabet::Alphabet;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::Deletion;
use crate::seq::indel::InDel;
use eyre::Report;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_primitives::Seq;
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqInfo {
  pub gaps: Vec<(usize, usize)>,
  pub unknown: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>,
  pub variable_indel: BTreeMap<(usize, usize), Deletion>,
  pub sequence: Seq,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DenseNodePartition {
  pub seq: DenseSeqInfo,
  pub profile: DenseSeqDistribution,
}

impl DenseNodePartition {
  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]);

    Ok(Self {
      seq: DenseSeqInfo {
        gaps,
        unknown,
        non_char,
        variable_indel: BTreeMap::new(),
        sequence: seq.to_owned(),
      },
      profile: DenseSeqDistribution::default(),
    })
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgePartition {
  pub indels: Vec<InDel>,
  pub transmission: Option<Vec<(usize, usize)>>,
  pub msg_to_child: DenseSeqDistribution,
  pub msg_to_parent: DenseSeqDistribution,
  pub msg_from_child: DenseSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqDistribution {
  pub dis: Array2<f64>,

  /// Total log likelihood
  pub log_lh: f64,
}

impl DenseSeqDistribution {
  pub fn new(dis: Array2<f64>, log_lh: f64) -> Self {
    Self { dis, log_lh }
  }
}
