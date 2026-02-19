use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use treetime_primitives::AsciiChar;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct Composition {
  counts: BTreeMap<AsciiChar, usize>,
  gap: AsciiChar,
}

impl Composition {
  /// Initialize counters with zeros, given an alphabet
  pub fn new<I>(alphabet_chars: I, gap: AsciiChar) -> Self
  where
    I: IntoIterator<Item = AsciiChar>,
  {
    let counts = alphabet_chars.into_iter().map(|c| (c, 0)).collect();
    Self { counts, gap }
  }

  /// Construct a `Composition` directly from a precomputed map and a gap character.
  pub fn from_counts<I: IntoIterator<Item = (AsciiChar, usize)>>(counts: I, gap: AsciiChar) -> Self {
    Self {
      counts: counts.into_iter().collect(),
      gap,
    }
  }

  pub fn get(&self, c: AsciiChar) -> Option<usize> {
    self.counts.get(&c).copied()
  }

  pub fn counts(&self) -> &BTreeMap<AsciiChar, usize> {
    &self.counts
  }

  /// Initialize counters to the composition of a given sequence
  pub fn with_sequence<SI, AI>(sequence: SI, alphabet_chars: AI, gap: AsciiChar) -> Self
  where
    SI: IntoIterator<Item = AsciiChar>,
    AI: IntoIterator<Item = AsciiChar>,
  {
    let mut this = Self::new(alphabet_chars, gap);
    this.add_sequence(sequence);
    this
  }

  /// Add composition of a given sequence to the counts
  pub fn add_sequence<I>(&mut self, sequence: I)
  where
    I: IntoIterator<Item = AsciiChar>,
  {
    let counts = sequence.into_iter().counts();
    self.counts.extend(counts);
  }

  /// Reflect sequence mutation in the composition counts
  pub fn add_sub(&mut self, sub: &Sub) {
    self.adjust_count(sub.reff(), -1);
    self.adjust_count(sub.qry(), 1);
  }

  /// Reflect sequence indel in the composition counts
  pub fn add_indel(&mut self, indel: &InDel) {
    let adjust_by = if indel.deletion { -1 } else { 1 };
    for nuc in &indel.seq {
      self.adjust_count(*nuc, adjust_by);
      self.adjust_count(self.gap, -adjust_by);
    }
  }

  pub fn adjust_count(&mut self, nuc: AsciiChar, change: isize) {
    let count = self.counts.entry(nuc).or_default();
    *count = count.saturating_add_signed(change);
  }
}
