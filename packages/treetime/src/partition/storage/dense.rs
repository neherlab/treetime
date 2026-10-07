use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::shared::update::MarginalNodeState;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::InDel;
use ndarray::Array2;
use std::collections::BTreeSet;
use treetime_primitives::{LogLh, Seq};
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug)]
pub struct DenseNodeState {
  pub seq: DenseSeqInfo,
  pub profile: DenseSeqDistribution,
}

#[derive(Clone, Debug)]
pub(crate) struct DenseLeafObs {
  sequence: Seq,
  gaps: Vec<(usize, usize)>,
  unknown: Vec<(usize, usize)>,
  non_char: Vec<(usize, usize)>,
}

impl DenseLeafObs {
  pub(crate) fn new(seq: &Seq, alphabet: &Alphabet) -> Self {
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]);
    Self {
      sequence: seq.to_owned(),
      gaps,
      unknown,
      non_char,
    }
  }

  pub(crate) fn sequence(&self) -> &Seq {
    &self.sequence
  }

  pub(crate) fn seq_info(&self) -> DenseSeqInfo {
    DenseSeqInfo {
      gaps: self.gaps.clone(),
      unknown: self.unknown.clone(),
      non_char: self.non_char.clone(),
      variable_indel: BTreeSet::new(),
      sequence: self.sequence.clone(),
    }
  }
}

impl MarginalNodeState for DenseNodeState {
  fn log_lh(&self) -> LogLh {
    self.profile.log_lh
  }
}

#[derive(Clone, Default, Debug)]
pub struct DenseSeqInfo {
  pub(crate) gaps: Vec<(usize, usize)>,
  pub(crate) unknown: Vec<(usize, usize)>,
  pub(crate) non_char: Vec<(usize, usize)>,
  pub(crate) variable_indel: BTreeSet<(usize, usize)>,
  pub(crate) sequence: Seq,
}

#[derive(Clone, Default, Debug)]
pub struct DenseEdgeBackward {
  pub(crate) msg_to_parent: DenseSeqDistribution,
  pub(crate) msg_from_child: DenseSeqDistribution,
}

#[derive(Clone, Default, Debug)]
pub struct DenseEdgeForward {
  pub(crate) msg_to_child: DenseSeqDistribution,
}

#[derive(Clone, Default, Debug)]
pub struct DenseEdgeEstimate {
  pub(crate) indels: Vec<InDel>,
}

#[derive(Clone, Default, Debug)]
pub struct DenseSeqDistribution {
  pub(crate) dis: Array2<f64>,

  pub(crate) log_lh: LogLh,
}

impl DenseSeqDistribution {
  pub fn new(dis: Array2<f64>, log_lh: LogLh) -> Self {
    Self { dis, log_lh }
  }
}
