use crate::alphabet::alphabet::Alphabet;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::InDel;
use eyre::Report;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use treetime_primitives::Seq;

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqInfo {
  pub gaps: Vec<(usize, usize)>,
  pub sequence: Seq,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DenseNodePartition {
  pub seq: DenseSeqInfo,
  pub profile: DenseSeqDis,
}

impl DenseNodePartition {
  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    let gaps = find_letter_ranges(seq, alphabet.gap());

    Ok(Self {
      seq: DenseSeqInfo {
        gaps,
        sequence: seq.to_owned(), // TODO(perf): try to avoid cloning
      },
      profile: DenseSeqDis::default(),
    })
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgePartition {
  pub indels: Vec<InDel>,
  pub transmission: Option<Vec<(usize, usize)>>,
  pub msg_to_child: DenseSeqDis,
  pub msg_to_parent: DenseSeqDis,
  pub msg_from_child: DenseSeqDis,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqDis {
  pub dis: Array2<f64>,

  /// Total log likelihood
  pub log_lh: f64,
}

impl DenseSeqDis {
  pub fn new(dis: Array2<f64>) -> Self {
    Self { dis, log_lh: 0.0 }
  }
}
