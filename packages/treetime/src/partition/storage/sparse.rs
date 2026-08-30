use crate::alphabet::alphabet::Alphabet;
use crate::seq::composition::Composition;
use crate::seq::find_char_ranges::find_letter_ranges;
use crate::seq::indel::{InDel, compose_indels, sort_indels};
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::mem;
use treetime_primitives::AlphabetLike;
use treetime_primitives::{AsciiChar, LogLh, Seq, StateSet, seq};
use treetime_utils::interval::range_union::range_union;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseNodePartition {
  pub seq: SparseSeqInfo,
  pub profile: SparseSeqDistribution,
}

impl SparseNodePartition {
  /// Create an empty placeholder node partition for edge split operations.
  /// The sequence and profile will be computed during the marginal update pass.
  pub fn empty(alphabet: &Alphabet) -> Self {
    Self {
      seq: SparseSeqInfo {
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
      profile: SparseSeqDistribution::default(),
    }
  }

  pub fn new(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    let variable = seq
      .iter()
      .enumerate()
      .filter(|&(_, c)| alphabet.is_ambiguous(*c))
      .map(|(pos, &c)| (pos, alphabet.char_to_set(c)))
      .collect();

    let seq_dis = FitchSeqDistribution {
      variable,
      variable_indel: BTreeSet::new(),
      chosen_state: btreemap! {},
    };

    let unknown = find_letter_ranges(seq, alphabet.unknown());
    let gaps = find_letter_ranges(seq, alphabet.gap());
    let non_char = range_union(&[unknown.clone(), gaps.clone()]); // TODO(perf): avoid cloning

    Ok(Self {
      seq: SparseSeqInfo {
        unknown,
        gaps,
        non_char,
        composition: Composition::with_seq(seq, alphabet.chars(), alphabet.gap()),
        sequence: seq.to_owned(), // TODO(perf): try to avoid cloning
        fitch: seq_dis,
      },
      profile: SparseSeqDistribution {
        variable: btreemap! {},
        variable_indel: BTreeSet::new(),
        fixed: btreemap! {},
        fixed_counts: Composition::new(alphabet.chars(), alphabet.gap()),
        log_lh: LogLh::ZERO,
      },
    })
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqInfo {
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub composition: Composition,      // count of all characters in the region that is not `non_char`
  pub sequence: Seq,
  pub fitch: FitchSeqDistribution,
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
#[allow(clippy::partial_pub_fields)]
pub struct SparseEdgePartition {
  subs_fitch: Vec<Sub>,
  subs_ml: Option<Vec<Sub>>,
  pub indels: Vec<InDel>,
  pub msg_to_parent: SparseSeqDistribution,
  pub msg_to_child: SparseSeqDistribution,
  pub msg_from_child: SparseSeqDistribution,

  /// The parent posterior evolved across this branch to the child (the marginal down-message).
  ///
  /// Populated by the forward pass only for edges whose child is a leaf, where tip imputation needs it
  /// at reconstruction time. `msg_to_child` stores the pre-propagation cavity message; the forward pass
  /// already computes this propagated form for the profile update but otherwise discards it. At a tip
  /// position with no observed state the leaf likelihood is uniform, so this down-message is the leaf's
  /// marginal posterior there, matching v0's per-leaf marginal profile.
  #[serde(default)]
  pub msg_from_parent: SparseSeqDistribution,

  pub transmission: Option<Vec<(usize, usize)>>,
}

impl SparseEdgePartition {
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
    self.subs_ml = None;
  }

  pub fn extend_fitch_subs(&mut self, subs: impl IntoIterator<Item = Sub>) {
    self.subs_fitch.extend(subs);
    self.subs_ml = None;
  }

  pub fn invert_fitch_subs(&mut self) {
    for sub in &mut self.subs_fitch {
      sub.invert();
    }
    self.subs_ml = None;
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

  pub fn ml_subs(&self) -> Option<&[Sub]> {
    self.subs_ml.as_deref()
  }

  pub fn set_ml_subs(&mut self, subs: Vec<Sub>) {
    self.subs_ml = Some(subs);
  }

  pub fn clear_ml_subs(&mut self) {
    self.subs_ml = None;
  }

  pub fn invert_for_reroot(&mut self) {
    self.invert_fitch_subs();
    self.clear_ml_subs();
    for indel in &mut self.indels {
      indel.invert();
    }
    mem::swap(&mut self.msg_to_parent, &mut self.msg_to_child);
    self.msg_from_child = SparseSeqDistribution::default();
    self.msg_from_parent = SparseSeqDistribution::default();
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseSeqDistribution {
  /// probability vector for each variable position collecting information from children
  pub variable: BTreeMap<usize, VarPos>,

  pub variable_indel: BTreeSet<(usize, usize)>,

  /// probability vector for the state of fixed positions based on information from children
  pub fixed: BTreeMap<AsciiChar, Array1<f64>>,

  pub fixed_counts: Composition,

  /// Total log likelihood
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
pub struct VarPos {
  pub dis: Array1<f64>, // array of floats of size 'alphabet'
  pub state: AsciiChar, // exact reference state for this sparse position
}

impl VarPos {
  pub fn new(dis: Array1<f64>, state: AsciiChar) -> Self {
    Self { dis, state }
  }
}
