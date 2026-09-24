use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::shared::update::MarginalNodeState;
use crate::seq::composition::Composition;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::{InDel, compose_indels, sort_indels};
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use treetime_primitives::AlphabetLike;
use treetime_primitives::{AsciiChar, LogLh, Seq, StateSet, seq};
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseNodeObs {
  pub(crate) unknown: Vec<(usize, usize)>,
  pub(crate) gaps: Vec<(usize, usize)>,
  pub(crate) non_char: Vec<(usize, usize)>,
  pub(crate) composition: Composition,
  pub(crate) fitch: FitchSeqDistribution,
}

impl SparseNodeObs {
  pub(crate) fn empty(alphabet: &Alphabet) -> Self {
    Self {
      unknown: vec![],
      gaps: vec![],
      non_char: vec![],
      composition: Composition::new(alphabet.chars(), alphabet.gap()),
      fitch: FitchSeqDistribution {
        variable: btreemap! {},
        variable_indel: BTreeSet::new(),
        chosen_state: btreemap! {},
      },
    }
  }

  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Self {
    let variable = seq
      .iter()
      .enumerate()
      .filter(|&(_, c)| alphabet.is_ambiguous(*c))
      .map(|(pos, &c)| (pos, alphabet.char_to_set(c)))
      .collect();

    let fitch = FitchSeqDistribution {
      variable,
      variable_indel: BTreeSet::new(),
      chosen_state: btreemap! {},
    };

    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]);

    Self {
      unknown,
      gaps,
      non_char,
      composition: Composition::with_seq(seq, alphabet.chars(), alphabet.gap()),
      fitch,
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseNodeState {
  pub(crate) sequence: Seq,
  pub(crate) profile: SparseSeqDistribution,

  #[serde(default)]
  pub(crate) emitted: Option<Seq>,
}

impl SparseNodeState {
  pub(crate) fn empty() -> Self {
    Self {
      sequence: seq![],
      profile: SparseSeqDistribution::default(),
      emitted: None,
    }
  }

  pub fn leaf(seq: &Seq) -> Self {
    Self {
      sequence: seq.to_owned(),
      profile: SparseSeqDistribution::default(),
      emitted: None,
    }
  }
}

impl MarginalNodeState for SparseNodeState {
  fn log_lh(&self) -> LogLh {
    self.profile.log_lh
  }

  fn set_log_lh(&mut self, log_lh: LogLh) {
    self.profile.log_lh = log_lh;
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
#[expect(
  clippy::partial_pub_fields,
  reason = "private fields hold derived state that only the constructor keeps consistent"
)]
pub struct SparseEdgeObs {
  subs_fitch: Vec<Sub>,
  pub indels: Vec<InDel>,
  pub(crate) transmission: Option<Vec<(usize, usize)>>,
}

impl SparseEdgeObs {
  pub fn with_fitch_subs(subs: Vec<Sub>) -> Self {
    Self {
      subs_fitch: subs,
      ..Default::default()
    }
  }

  pub(crate) fn fitch_subs(&self) -> &[Sub] {
    &self.subs_fitch
  }

  pub(crate) fn set_fitch_subs(&mut self, subs: Vec<Sub>) {
    self.subs_fitch = subs;
  }

  pub(crate) fn extend_fitch_subs(&mut self, subs: impl IntoIterator<Item = Sub>) {
    self.subs_fitch.extend(subs);
  }

  pub(crate) fn invert_fitch_subs(&mut self) {
    for sub in &mut self.subs_fitch {
      sub.invert();
    }
  }

  pub(crate) fn chain_fitch_subs(&self, suffix: &[Sub]) -> Result<Vec<Sub>, Report> {
    compose_substitutions(&self.subs_fitch, suffix)
  }

  pub(crate) fn chain_fitch_indels(&self, child_indels: &[InDel]) -> Vec<InDel> {
    let mut parent = self.indels.clone();
    let mut child = child_indels.to_vec();
    sort_indels(&mut parent);
    sort_indels(&mut child);
    compose_indels(&parent, &child)
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgeBackward {
  pub(crate) msg_to_parent: SparseSeqDistribution,
  pub(crate) msg_from_child: SparseSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgeForward {
  pub(crate) msg_to_child: SparseSeqDistribution,

  #[serde(default)]
  pub(crate) msg_from_parent: SparseSeqDistribution,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqDistribution {
  pub(crate) variable: BTreeMap<usize, VarPos>,

  pub(crate) variable_indel: BTreeSet<(usize, usize)>,

  pub(crate) fixed: BTreeMap<AsciiChar, Array1<f64>>,

  pub(crate) fixed_counts: Composition,

  pub(crate) log_lh: LogLh,
}

impl Default for SparseSeqDistribution {
  fn default() -> Self {
    Self {
      variable: btreemap! {},
      variable_indel: BTreeSet::new(),
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: LogLh::ZERO,
    }
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitchNodeData {
  pub(crate) seq: FitchSeqInfo,
}

impl FitchNodeData {
  pub(crate) fn empty(alphabet: &Alphabet) -> Self {
    Self {
      seq: FitchSeqInfo {
        unknown: vec![],
        gaps: vec![],
        non_char: vec![],
        composition: Composition::new(alphabet.chars(), alphabet.gap()),
        sequence: seq![],
        fitch: FitchSeqDistribution {
          variable: btreemap! {},
          variable_indel: BTreeSet::new(),
          chosen_state: btreemap! {},
        },
      },
    }
  }

  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    let variable = seq
      .iter()
      .enumerate()
      .filter(|&(_, c)| alphabet.is_ambiguous(*c))
      .map(|(pos, &c)| (pos, alphabet.char_to_set(c)))
      .collect();

    let fitch = FitchSeqDistribution {
      variable,
      variable_indel: BTreeSet::new(),
      chosen_state: btreemap! {},
    };

    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]);

    Ok(Self {
      seq: FitchSeqInfo {
        unknown,
        gaps,
        non_char,
        composition: Composition::with_seq(seq, alphabet.chars(), alphabet.gap()),
        sequence: seq.to_owned(),
        fitch,
      },
    })
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitchSeqInfo {
  pub(crate) unknown: Vec<(usize, usize)>,
  pub(crate) gaps: Vec<(usize, usize)>,
  pub(crate) non_char: Vec<(usize, usize)>,
  pub(crate) composition: Composition,
  pub(crate) sequence: Seq,
  pub(crate) fitch: FitchSeqDistribution,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitchSeqDistribution {
  pub(crate) variable: BTreeMap<usize, StateSet>,

  pub(crate) variable_indel: BTreeSet<(usize, usize)>,

  pub(crate) chosen_state: BTreeMap<usize, AsciiChar>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VarPos {
  pub(crate) dis: Array1<f64>,
  pub(crate) state: AsciiChar,
}

impl VarPos {}
