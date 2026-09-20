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
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>,
  pub composition: Composition,
  pub fitch: FitchSeqDistribution,
}

impl SparseNodeObs {
  pub fn empty(alphabet: &Alphabet) -> Self {
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
  pub sequence: Seq,
  pub profile: SparseSeqDistribution,

  #[serde(default)]
  pub emitted: Option<Seq>,
}

impl SparseNodeState {
  pub fn empty() -> Self {
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
#[allow(clippy::partial_pub_fields)]
pub struct SparseEdgeObs {
  subs_fitch: Vec<Sub>,
  pub indels: Vec<InDel>,
  pub transmission: Option<Vec<(usize, usize)>>,
}

impl SparseEdgeObs {
  pub fn with_fitch_subs(subs: Vec<Sub>) -> Self {
    Self {
      subs_fitch: subs,
      ..Default::default()
    }
  }

  pub fn with_fitch_subs_and_indels(subs: Vec<Sub>, indels: Vec<InDel>) -> Self {
    Self {
      subs_fitch: subs,
      indels,
      ..Default::default()
    }
  }

  pub fn fitch_subs(&self) -> &[Sub] {
    &self.subs_fitch
  }

  pub fn set_fitch_subs(&mut self, subs: Vec<Sub>) {
    self.subs_fitch = subs;
  }

  pub fn extend_fitch_subs(&mut self, subs: impl IntoIterator<Item = Sub>) {
    self.subs_fitch.extend(subs);
  }

  pub fn invert_fitch_subs(&mut self) {
    for sub in &mut self.subs_fitch {
      sub.invert();
    }
  }

  pub fn chain_fitch_subs(&self, suffix: &[Sub]) -> Result<Vec<Sub>, Report> {
    compose_substitutions(&self.subs_fitch, suffix)
  }

  pub fn chain_fitch_indels(&self, child_indels: &[InDel]) -> Vec<InDel> {
    let mut parent = self.indels.clone();
    let mut child = child_indels.to_vec();
    sort_indels(&mut parent);
    sort_indels(&mut child);
    compose_indels(&parent, &child)
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgeBackward {
  pub msg_to_parent: SparseSeqDistribution,
  pub msg_from_child: SparseSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgeForward {
  pub msg_to_child: SparseSeqDistribution,

  #[serde(default)]
  pub msg_from_parent: SparseSeqDistribution,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqDistribution {
  pub variable: BTreeMap<usize, VarPos>,

  pub variable_indel: BTreeSet<(usize, usize)>,

  pub fixed: BTreeMap<AsciiChar, Array1<f64>>,

  pub fixed_counts: Composition,

  pub log_lh: LogLh,
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
pub struct FitchSeqDistribution {
  pub variable: BTreeMap<usize, StateSet>,

  pub variable_indel: BTreeSet<(usize, usize)>,

  pub chosen_state: BTreeMap<usize, AsciiChar>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FitchNodeData {
  pub seq: FitchSeqInfo,
}

impl FitchNodeData {
  pub fn empty(alphabet: &Alphabet) -> Self {
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
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>,
  pub composition: Composition,
  pub sequence: Seq,
  pub fitch: FitchSeqDistribution,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct VarPos {
  pub dis: Array1<f64>,
  pub state: AsciiChar,
}

impl VarPos {
  pub fn new(dis: Array1<f64>, state: AsciiChar) -> Self {
    Self { dis, state }
  }
}
