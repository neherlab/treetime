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
  pub fn new<I, IC, GC>(alphabet_chars: I, gap: GC) -> Self
  where
    I: IntoIterator<Item = IC>,
    IC: Into<AsciiChar>,
    GC: Into<AsciiChar>,
  {
    let counts = alphabet_chars.into_iter().map(|c| (c.into(), 0)).collect();
    let gap = gap.into();
    Self { counts, gap }
  }

  /// Construct a `Composition` directly from a precomputed map and a gap character.
  pub fn from<C: Into<AsciiChar>, I: IntoIterator<Item = (C, usize)>>(counts: I, gap: C) -> Self {
    Self {
      counts: counts.into_iter().map(|(c, count)| (c.into(), count)).collect(),
      gap: gap.into(),
    }
  }

  pub fn get<C: Into<AsciiChar>>(&self, c: C) -> Option<usize> {
    self.counts.get(&c.into()).copied()
  }

  pub fn counts(&self) -> &BTreeMap<AsciiChar, usize> {
    &self.counts
  }

  /// Initialize counters to the composition of a given sequence
  pub fn with_sequence<SC, SI, AI, AC, GC>(sequence: SI, alphabet_chars: AI, gap: GC) -> Self
  where
    SI: IntoIterator<Item = SC>,
    SC: Into<AsciiChar>,
    AI: IntoIterator<Item = AC>,
    AC: Into<AsciiChar>,
    GC: Into<AsciiChar>,
  {
    let mut this = Self::new(alphabet_chars, gap);
    this.add_sequence(sequence.into_iter().map(Into::into));
    this
  }

  /// Add composition of a given sequence to the counts
  pub fn add_sequence<C, I>(&mut self, sequence: I)
  where
    I: IntoIterator<Item = C>,
    C: Into<AsciiChar>,
  {
    let counts = sequence.into_iter().map(Into::into).counts();
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

  pub fn adjust_count<C: Into<AsciiChar>>(&mut self, nuc: C, change: isize) {
    let count = self.counts.entry(nuc.into()).or_default();
    *count = count.saturating_add_signed(change);
  }
}
