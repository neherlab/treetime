use crate::alphabet::alphabet::Alphabet;
use crate::seq::find_char_ranges::find_letter_ranges_by;
use itertools::Itertools;
use treetime_primitives::{AsciiChar, Seq};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Query {
  diffs: Vec<(usize, AsciiChar)>,
  missing: Vec<(usize, usize)>,
}

impl Query {
  pub fn from_seq(seq: &Seq, reference: &Seq, alphabet: &Alphabet) -> Self {
    let diffs = seq
      .as_slice()
      .iter()
      .zip(reference.as_slice())
      .enumerate()
      .filter(|&(_, (&state, &reff))| alphabet.is_canonical(state) && state != reff)
      .map(|(pos, (&state, _))| (pos, state))
      .collect_vec();
    let missing = find_letter_ranges_by(seq, |state| !alphabet.is_canonical(state));
    Self { diffs, missing }
  }

  pub(crate) fn from_diffs(diffs: Vec<(usize, AsciiChar)>) -> Self {
    Self { diffs, missing: vec![] }
  }

  pub(crate) fn diffs(&self) -> &[(usize, AsciiChar)] {
    &self.diffs
  }

  pub(crate) fn is_missing(&self, pos: usize) -> bool {
    let index = self.missing.partition_point(|&(_, end)| end <= pos);
    self.missing.get(index).is_some_and(|&(start, _)| start <= pos)
  }

  pub(crate) fn state(&self, pos: usize, reference: &Seq) -> Option<AsciiChar> {
    if let Ok(index) = self.diffs.binary_search_by_key(&pos, |&(p, _)| p) {
      return Some(self.diffs[index].1);
    }
    (!self.is_missing(pos)).then(|| reference[pos])
  }

  pub(crate) fn n_missing_in(&self, positions: &[usize]) -> usize {
    self
      .missing
      .iter()
      .map(|&(start, end)| positions.partition_point(|&pos| pos < end) - positions.partition_point(|&pos| pos < start))
      .sum()
  }

  pub fn n_missing(&self) -> usize {
    self.missing.iter().map(|(start, end)| end - start).sum()
  }

  pub(crate) fn cost_against(&self, node_state: &[(usize, AsciiChar)], reference: &Seq) -> usize {
    let mut cost = 0;
    for entry in self.diffs.iter().merge_join_by(node_state, |a, b| a.0.cmp(&b.0)) {
      let differs = match entry {
        itertools::EitherOrBoth::Both(&(_, query), &(_, node)) => query != node,
        itertools::EitherOrBoth::Left(&(pos, query)) => query != reference[pos],
        itertools::EitherOrBoth::Right(&(pos, node)) => !self.is_missing(pos) && reference[pos] != node,
      };
      cost += usize::from(differs);
    }
    cost
  }
}
