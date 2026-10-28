use crate::alphabet::alphabet::Alphabet;
use crate::partition::marginal::shared::update::MarginalNodeState;
use crate::partition::storage::var_pos_map::VarPosMap;
use crate::seq::composition::Composition;
use crate::seq::indel::InDel;
use crate::seq::mutation::{Sub, compose_substitutions};
use eyre::Report;
use maplit::btreemap;
use smallvec::SmallVec;
use std::collections::{BTreeMap, BTreeSet};
use std::ops::Index;
use std::sync::Arc;
use treetime_primitives::{AsciiChar, AsciiTable, LogLh, Seq, StateSet, seq};
use treetime_utils::interval::range_union::range_union;

pub(crate) const INLINE_STATES: usize = 4;

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
      composition: Composition::new(alphabet),
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

#[derive(Clone, Default, Debug, PartialEq, Eq)]
#[expect(
  clippy::partial_pub_fields,
  reason = "private fields hold derived state that only the constructor keeps consistent"
)]
pub struct SparseEdgeObs {
  subs_fitch: Vec<Sub>,
  pub indels: Vec<InDel>,
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

#[derive(Clone, Debug, Default)]
pub struct SparseSeqDistribution {
  pub(crate) variable: VarPosMap,

  pub(crate) fixed: FixedProfiles,

  pub(crate) log_lh: LogLh,
}

#[derive(Clone, Debug, Default)]
pub struct FixedProfiles {
  entries: SmallVec<[FixedProfile; INLINE_STATES]>,
}

impl FixedProfiles {
  pub(crate) fn get(&self, state: AsciiChar) -> Option<&FixedProfile> {
    let index = self.entries.partition_point(|entry| entry.state < state);
    self.entries.get(index).filter(|entry| entry.state == state)
  }

  pub(crate) fn insert(&mut self, state: AsciiChar, count: usize, dis: SiteProfile) {
    let entry = FixedProfile { state, count, dis };
    if self.entries.last().is_none_or(|last| last.state < state) {
      self.entries.push(entry);
      return;
    }
    match self.entries.binary_search_by_key(&state, |entry| entry.state) {
      Ok(index) => self.entries[index] = entry,
      Err(index) => self.entries.insert(index, entry),
    }
  }

  pub(crate) fn iter(&self) -> impl Iterator<Item = &FixedProfile> {
    self.entries.iter()
  }

  pub(crate) fn len(&self) -> usize {
    self.entries.len()
  }

  pub(crate) fn map_dis(&self, transform: impl Fn(&SiteProfile) -> SiteProfile) -> Self {
    Self {
      entries: self
        .entries
        .iter()
        .map(|entry| FixedProfile {
          state: entry.state,
          count: entry.count,
          dis: transform(&entry.dis),
        })
        .collect(),
    }
  }
}

impl Index<AsciiChar> for FixedProfiles {
  type Output = FixedProfile;

  #[expect(
    clippy::expect_used,
    reason = "a state without a fixed profile is a programming error, as an absent key is for a map index"
  )]
  fn index(&self, state: AsciiChar) -> &FixedProfile {
    self.get(state).expect("state has no fixed profile")
  }
}

#[derive(Clone, Debug)]
pub struct FixedProfile {
  pub(crate) state: AsciiChar,
  pub(crate) count: usize,
  pub(crate) dis: SiteProfile,
}

pub(crate) type SiteProfile = SmallVec<[f64; INLINE_STATES]>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FitchNodeData {
  pub(crate) seq: FitchSeqInfo,
}

impl FitchNodeData {
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FitchSeqInfo {
  pub(crate) unknown: Vec<(usize, usize)>,
  pub(crate) gaps: Vec<(usize, usize)>,
  pub(crate) non_char: Vec<(usize, usize)>,
  pub(crate) composition: Composition,
  pub(crate) sequence: Seq,
  pub(crate) fitch: FitchSeqDistribution,
}

#[derive(Clone, Copy, Debug)]
pub struct VarPos<'a> {
  pub(crate) dis: &'a [f64],
  pub(crate) state: AsciiChar,
}

fn observe_seq(seq: &Seq, alphabet: &Alphabet) -> SeqObservation {
  let flags = char_flags(alphabet);
  let mut variable = Vec::new();
  let mut unknown = RangeRuns::default();
  let mut gaps = RangeRuns::default();
  let mut previous = 0_u8;
  for (pos, &c) in seq.iter().enumerate() {
    let current = flags[c];
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
    composition: Composition::with_seq(seq, alphabet),
    fitch: FitchSeqDistribution {
      variable: variable.into_iter().collect(),
      variable_indel: BTreeSet::new(),
      chosen_state: btreemap! {},
    },
  }
}

struct SeqObservation {
  unknown: Vec<(usize, usize)>,
  gaps: Vec<(usize, usize)>,
  non_char: Vec<(usize, usize)>,
  composition: Composition,
  fitch: FitchSeqDistribution,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FitchSeqDistribution {
  pub(crate) variable: BTreeMap<usize, StateSet>,

  pub(crate) variable_indel: BTreeSet<(usize, usize)>,

  pub(crate) chosen_state: BTreeMap<usize, AsciiChar>,
}

fn char_flags(alphabet: &Alphabet) -> AsciiTable<u8> {
  AsciiTable::from_fn(|c| {
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
