use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::shared::update::MarginalNodeState;
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use crate::seq::indel::{InDel, compose_indels, sort_indels};
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_primitives::AlphabetLike;
use treetime_primitives::{AsciiChar, LogLh, Seq, StateSet, seq};
use treetime_utils::interval::range_union::range_union;

const AMBIGUOUS: u8 = 1;
const UNKNOWN: u8 = 2;
const GAP: u8 = 4;

#[derive(Clone, Debug)]
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
    let SeqObservation {
      unknown,
      gaps,
      non_char,
      composition,
      fitch,
    } = observe_seq(seq, alphabet);
    Self {
      unknown,
      gaps,
      non_char,
      composition,
      fitch,
    }
  }
}

#[derive(Clone, Debug)]
pub struct SparseNodeState {
  pub(crate) sequence: Arc<Seq>,
  pub(crate) profile: SparseSeqDistribution,
}

impl SparseNodeState {
  pub(crate) fn empty() -> Self {
    Self {
      sequence: Arc::new(seq![]),
      profile: SparseSeqDistribution::default(),
    }
  }

  pub fn leaf(seq: &Seq) -> Self {
    Self {
      sequence: Arc::new(seq.to_owned()),
      profile: SparseSeqDistribution::default(),
    }
  }
}

impl MarginalNodeState for SparseNodeState {
  fn log_lh(&self) -> LogLh {
    self.profile.log_lh
  }
}

#[derive(Clone, Default, Debug)]
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

#[derive(Clone, Default, Debug)]
pub struct SparseEdgeBackward {
  pub(crate) msg_to_parent: SparseSeqDistribution,
  pub(crate) msg_from_child: SparseSeqDistribution,
}

#[derive(Clone, Default, Debug)]
pub struct SparseEdgeForward {
  pub(crate) msg_to_child: SparseSeqDistribution,

  pub(crate) msg_from_parent: SparseSeqDistribution,
}

#[derive(Clone, Debug)]
pub struct SparseSeqDistribution {
  pub(crate) variable: VarPosMap,

  pub(crate) fixed: BTreeMap<AsciiChar, Array1<f64>>,

  pub(crate) fixed_counts: Composition,

  pub(crate) log_lh: LogLh,
}

impl Default for SparseSeqDistribution {
  fn default() -> Self {
    Self {
      variable: VarPosMap::default(),
      fixed: btreemap! {},
      fixed_counts: Composition::new(std::iter::empty::<AsciiChar>(), AsciiChar::from_byte_unchecked(b'-')),
      log_lh: LogLh::ZERO,
    }
  }
}

#[derive(Clone, Debug)]
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

  pub fn new(seq: Seq, alphabet: &Alphabet) -> Self {
    let SeqObservation {
      unknown,
      gaps,
      non_char,
      composition,
      fitch,
    } = observe_seq(&seq, alphabet);
    Self {
      seq: FitchSeqInfo {
        unknown,
        gaps,
        non_char,
        composition,
        sequence: seq,
        fitch,
      },
    }
  }
}

#[derive(Clone, Debug)]
pub struct FitchSeqInfo {
  pub(crate) unknown: Vec<(usize, usize)>,
  pub(crate) gaps: Vec<(usize, usize)>,
  pub(crate) non_char: Vec<(usize, usize)>,
  pub(crate) composition: Composition,
  pub(crate) sequence: Seq,
  pub(crate) fitch: FitchSeqDistribution,
}

#[derive(Clone, Debug)]
pub struct FitchSeqDistribution {
  pub(crate) variable: BTreeMap<usize, StateSet>,

  pub(crate) variable_indel: BTreeSet<(usize, usize)>,

  pub(crate) chosen_state: BTreeMap<usize, AsciiChar>,
}

#[derive(Clone, Debug)]
pub struct VarPos {
  pub(crate) dis: Array1<f64>,
  pub(crate) state: AsciiChar,
}

impl VarPos {}

struct SeqObservation {
  unknown: Vec<(usize, usize)>,
  gaps: Vec<(usize, usize)>,
  non_char: Vec<(usize, usize)>,
  composition: Composition,
  fitch: FitchSeqDistribution,
}

fn observe_seq(seq: &Seq, alphabet: &Alphabet) -> SeqObservation {
  let flags = char_flags(alphabet);
  let mut variable = Vec::new();
  let mut unknown = RangeRuns::default();
  let mut gaps = RangeRuns::default();
  let mut previous = 0_u8;
  for (pos, &c) in seq.iter().enumerate() {
    let current = flags[usize::from(c)];
    if current & AMBIGUOUS != 0 {
      variable.push((pos, alphabet.char_to_set(c)));
    }
    if current != previous {
      unknown.step(pos, current & UNKNOWN != 0);
      gaps.step(pos, current & GAP != 0);
      previous = current;
    }
  }
  let unknown = unknown.finish(seq.len());
  let gaps = gaps.finish(seq.len());
  let non_char = range_union(&[unknown.clone(), gaps.clone()]);
  SeqObservation {
    unknown,
    gaps,
    non_char,
    composition: Composition::with_seq(seq, alphabet.chars(), alphabet.gap()),
    fitch: FitchSeqDistribution {
      variable: variable.into_iter().collect(),
      variable_indel: BTreeSet::new(),
      chosen_state: btreemap! {},
    },
  }
}

#[allow(
  clippy::as_conversions,
  reason = "table indices below 128 convert exactly to ASCII bytes"
)]
fn char_flags(alphabet: &Alphabet) -> [u8; 128] {
  std::array::from_fn(|index| {
    let c = AsciiChar::from_byte_unchecked(index as u8);
    let mut flags = 0;
    if alphabet.is_ambiguous(c) {
      flags |= AMBIGUOUS;
    }
    if c == alphabet.unknown() {
      flags |= UNKNOWN;
    }
    if c == alphabet.gap() {
      flags |= GAP;
    }
    flags
  })
}

#[derive(Debug, Default)]
struct RangeRuns {
  ranges: Vec<(usize, usize)>,
  start: Option<usize>,
}

impl RangeRuns {
  fn step(&mut self, pos: usize, inside: bool) {
    match (inside, self.start) {
      (true, None) => self.start = Some(pos),
      (false, Some(start)) => {
        self.ranges.push((start, pos));
        self.start = None;
      },
      _ => {},
    }
  }

  fn finish(mut self, length: usize) -> Vec<(usize, usize)> {
    if let Some(start) = self.start {
      self.ranges.push((start, length));
    }
    self.ranges
  }
}
