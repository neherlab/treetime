use crate::alphabet::alphabet::Alphabet;
use crate::seq::find_char_ranges::find_letter_ranges_by;
use itertools::{EitherOrBoth, Itertools};
use treetime_primitives::{AsciiChar, Seq, StateSet};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Query {
  diffs: Vec<(usize, StateSet)>,
  missing: Vec<(usize, usize)>,
}

impl Query {
  pub fn from_seq(seq: &Seq, reference: &Seq, alphabet: &Alphabet) -> Self {
    let states = |state: AsciiChar| {
      (!alphabet.is_unknown(state) && !alphabet.is_gap(state))
        .then(|| alphabet.canonical_states(state).ok())
        .flatten()
        .filter(|states| !states.is_empty())
    };
    let diffs = seq
      .as_slice()
      .iter()
      .zip(reference.as_slice())
      .enumerate()
      .filter(|&(_, (&state, &reff))| state != reff)
      .filter_map(|(pos, (&state, _))| Some((pos, states(state)?)))
      .collect_vec();
    let missing = find_letter_ranges_by(seq, |state| states(state).is_none());
    Self { diffs, missing }
  }

  pub(crate) fn from_state_and_sets(state: &[(usize, AsciiChar)], sets: &[(usize, StateSet)], reference: &Seq) -> Self {
    let diffs = state
      .iter()
      .map(|&(pos, value)| (pos, StateSet::from_char(value)))
      .merge_join_by(sets, |a, b| a.0.cmp(&b.0))
      .map(|entry| match entry {
        EitherOrBoth::Left(single) => single,
        EitherOrBoth::Right(&set) | EitherOrBoth::Both(_, &set) => set,
      })
      .filter(|&(pos, set)| !(set.len() == 1 && set.contains(reference[pos])))
      .collect();
    Self { diffs, missing: vec![] }
  }

  pub(crate) fn is_missing(&self, pos: usize) -> bool {
    let index = self.missing.partition_point(|&(_, end)| end <= pos);
    self.missing.get(index).is_some_and(|&(start, _)| start <= pos)
  }

  pub(crate) fn states(&self, pos: usize, reference: &Seq) -> Option<StateSet> {
    if let Ok(index) = self.diffs.binary_search_by_key(&pos, |&(p, _)| p) {
      return Some(self.diffs[index].1);
    }
    (!self.is_missing(pos)).then(|| StateSet::from_char(reference[pos]))
  }

  pub(crate) fn entries(&self) -> &[(usize, StateSet)] {
    &self.diffs
  }

  pub fn n_missing(&self) -> usize {
    self.missing.iter().map(|(start, end)| end - start).sum()
  }

  pub(crate) fn root_cost(&self, reference: &Seq) -> usize {
    self
      .diffs
      .iter()
      .filter(|&&(pos, states)| !states.contains(reference[pos]))
      .count()
  }

  pub(crate) fn distance(&self, other: &Self, reference: &Seq) -> usize {
    let single = |pos: usize| StateSet::from_char(reference[pos]);
    self
      .diffs
      .iter()
      .merge_join_by(&other.diffs, |a, b| a.0.cmp(&b.0))
      .filter(|entry| {
        let (own, theirs) = match *entry {
          EitherOrBoth::Both(&(_, own), &(_, theirs)) => (own, theirs),
          EitherOrBoth::Left(&(pos, own)) => {
            if other.is_missing(pos) {
              return false;
            }
            (own, single(pos))
          },
          EitherOrBoth::Right(&(pos, theirs)) => {
            if self.is_missing(pos) {
              return false;
            }
            (single(pos), theirs)
          },
        };
        own.intersection(&theirs).is_empty()
      })
      .count()
  }

  pub(crate) fn cost_against(&self, node_state: &[(usize, AsciiChar)], reference: &Seq) -> usize {
    let mut cost = 0;
    for entry in self.diffs.iter().merge_join_by(node_state, |a, b| a.0.cmp(&b.0)) {
      let differs = match entry {
        EitherOrBoth::Both(&(_, states), &(_, node)) => !states.contains(node),
        EitherOrBoth::Left(&(pos, states)) => !states.contains(reference[pos]),
        EitherOrBoth::Right(&(pos, node)) => !self.is_missing(pos) && reference[pos] != node,
      };
      cost += usize::from(differs);
    }
    cost
  }
}
