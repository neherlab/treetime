use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::shared::update::MarginalNodeState;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::InDel;
use eyre::Report;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
use treetime_primitives::{LogLh, Seq};
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DenseNodeState {
  pub seq: DenseSeqInfo,
  pub profile: DenseSeqDistribution,
}

impl DenseNodeState {
  pub(crate) fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]);

    Ok(Self {
      seq: DenseSeqInfo {
        gaps,
        unknown,
        non_char,
        variable_indel: BTreeSet::new(),
        sequence: seq.to_owned(),
      },
      profile: DenseSeqDistribution::default(),
    })
  }

  pub(crate) fn empty() -> Self {
    Self {
      seq: DenseSeqInfo::default(),
      profile: DenseSeqDistribution::default(),
    }
  }
}

impl MarginalNodeState for DenseNodeState {
  fn log_lh(&self) -> LogLh {
    self.profile.log_lh
  }

  fn set_log_lh(&mut self, log_lh: LogLh) {
    self.profile.log_lh = log_lh;
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqInfo {
  pub(crate) gaps: Vec<(usize, usize)>,
  pub(crate) unknown: Vec<(usize, usize)>,
  pub(crate) non_char: Vec<(usize, usize)>,
  pub(crate) variable_indel: BTreeSet<(usize, usize)>,
  pub(crate) sequence: Seq,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgeBackward {
  pub(crate) msg_to_parent: DenseSeqDistribution,
  pub(crate) msg_from_child: DenseSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgeForward {
  pub(crate) msg_to_child: DenseSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgeEstimate {
  pub(crate) indels: Vec<InDel>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqDistribution {
  pub(crate) dis: Array2<f64>,

  pub(crate) log_lh: LogLh,
}

impl DenseSeqDistribution {
  pub fn new(dis: Array2<f64>, log_lh: LogLh) -> Self {
    Self { dis, log_lh }
  }
}
