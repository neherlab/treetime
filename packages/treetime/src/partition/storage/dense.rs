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

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqInfo {
  pub gaps: Vec<(usize, usize)>,
  pub unknown: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>,
  pub variable_indel: BTreeSet<(usize, usize)>,
  pub sequence: Seq,
}

/// The evolving per-node dense state: site information and the posterior profile.
///
/// This is a genuinely unified positional slot, not a stage split: `seq.sequence` is the observed
/// residue at a leaf and the reconstructed residue at an internal node, and `profile` is refined by the
/// backward and then the forward pass. Both passes borrow the durable inputs and return updated node
/// states; nothing here is a durable input the partition owns.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DenseNodeState {
  pub seq: DenseSeqInfo,
  pub profile: DenseSeqDistribution,
}

impl DenseNodeState {
  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
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

  pub fn empty() -> Self {
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

/// Backward-pass edge messages, produced by the backward pass and consumed by the forward pass and by
/// transition counting. Distinct owner from the forward messages and the final estimates.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgeBackward {
  pub msg_to_parent: DenseSeqDistribution,
  pub msg_from_child: DenseSeqDistribution,
}

/// Forward-pass edge message, produced by the forward pass and consumed by transition counting and the
/// branch-length optimizer. Distinct owner from the backward messages and the final estimates.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgeForward {
  pub msg_to_child: DenseSeqDistribution,
}

/// Final per-edge estimate produced by the forward pass: the indels placed on the branch.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseEdgeEstimate {
  pub indels: Vec<InDel>,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct DenseSeqDistribution {
  pub dis: Array2<f64>,

  /// Total log likelihood
  pub log_lh: LogLh,
}

impl DenseSeqDistribution {
  pub fn new(dis: Array2<f64>, log_lh: LogLh) -> Self {
    Self { dis, log_lh }
  }
}
