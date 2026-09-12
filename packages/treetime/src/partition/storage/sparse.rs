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
use treetime_primitives::AlphabetLike;
use treetime_primitives::{AsciiChar, LogLh, Seq, StateSet, seq};
use treetime_utils::interval::range_union::range_union;

/// Durable per-node observations produced by the Fitch pre-pass and consumed (never rewritten) by the
/// marginal passes and reconstruction: the ambiguity/gap ranges, the residue composition, and the
/// Fitch state sets. The partition owns these as immutable inputs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseNodeObs {
  pub unknown: Vec<(usize, usize)>,
  pub gaps: Vec<(usize, usize)>,
  pub non_char: Vec<(usize, usize)>, // any position that does not evolve according to the substitution model, i.e. gap or N
  pub composition: Composition,      // count of all characters in the region that is not `non_char`
  pub fitch: FitchSeqDistribution,
}

impl SparseNodeObs {
  /// Empty observations for a placeholder node created by an edge split.
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

/// The evolving per-node marginal state: the reconstructed (or observed, at a leaf) sequence, the
/// posterior profile, and the cached emitted output. A genuinely unified positional slot -- `sequence`
/// is observed at a leaf and parsimony-reconstructed at an internal node -- so it is kept together
/// rather than split across obs/result structs.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SparseNodeState {
  pub sequence: Seq,
  pub profile: SparseSeqDistribution,

  /// Cached output sequence for the two cases that cannot be rebuilt from `sequence` and `profile`.
  ///
  /// Every output path normally rebuilds the sequence on demand, so nothing needs storing. Two results
  /// cannot be rebuilt and are kept here so all paths return the same bytes: a random draw under
  /// `--sample-from-profile`, and an imputed leaf, whose filled-in states depend on
  /// `--impute-missing-data`. `None` in a default run.
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

  /// Seed the node state for a leaf from its observed sequence: the observed sequence is also the
  /// leaf's initial marginal sequence.
  pub fn leaf(seq: &Seq) -> Self {
    Self {
      sequence: seq.to_owned(),
      profile: SparseSeqDistribution::default(),
      emitted: None,
    }
  }
}

/// Durable per-edge observations produced by the Fitch pre-pass: the Fitch substitutions, the indels,
/// and the transmission mask. The partition owns these as immutable inputs to the marginal passes.
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

/// Backward-pass edge messages: distinct owner from the forward messages and the final estimates.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgeBackward {
  pub msg_to_parent: SparseSeqDistribution,
  pub msg_from_child: SparseSeqDistribution,
}

/// Forward-pass edge messages: the cavity down-message and, for leaf edges, the propagated parent
/// posterior used for tip imputation. Distinct owner from the backward messages and the estimates.
#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct SparseEdgeForward {
  pub msg_to_child: SparseSeqDistribution,

  /// The parent posterior evolved across this branch to the child (the marginal down-message).
  ///
  /// Populated by the forward pass only for edges whose child is a leaf, where tip imputation needs it
  /// at reconstruction time. `msg_to_child` stores the pre-propagation cavity message; the forward pass
  /// already computes this propagated form for the profile update but otherwise discards it. At a tip
  /// position with no observed state the leaf likelihood is uniform, so this down-message is the leaf's
  /// marginal posterior there, matching v0's per-leaf marginal profile.
  #[serde(default)]
  pub msg_from_parent: SparseSeqDistribution,
}

/// Final per-edge marginal estimate: the maximum-likelihood substitutions placed on the branch,
/// produced by the forward pass and held in the estimates map as a plain `Vec<Sub>`.
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

/// The Fitch parsimony pre-pass's own working per-node data. The parsimony passes build this in place;
/// the marginal handoff (`PartitionFitch::into_marginal_*`) splits it into the durable observations the
/// partition owns and the seed node state the marginal passes evolve.
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
  pub dis: Array1<f64>, // array of floats of size 'alphabet'
  pub state: AsciiChar, // exact reference state for this sparse position
}

impl VarPos {
  pub fn new(dis: Array1<f64>, state: AsciiChar) -> Self {
    Self { dis, state }
  }
}
