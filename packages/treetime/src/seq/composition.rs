use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use eyre::Report;
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

  pub fn with_seq(
    sequence: impl AsRef<[AsciiChar]>,
    alphabet_chars: impl IntoIterator<Item = AsciiChar>,
    gap: AsciiChar,
  ) -> Self {
    let mut this = Self::new(alphabet_chars, gap);
    this.add_seq(sequence);
    this
  }

  pub fn with_seq_str(
    sequence: &str,
    alphabet_chars: impl IntoIterator<Item = AsciiChar>,
    gap: AsciiChar,
  ) -> Result<Self, Report> {
    let seq = Self::str_to_chars(sequence)?;
    Ok(Self::with_seq(seq, alphabet_chars, gap))
  }

  pub fn add_seq(&mut self, sequence: impl AsRef<[AsciiChar]>) {
    let mut additions = [0; 128];
    for &c in sequence.as_ref() {
      additions[usize::from(c)] += 1;
    }
    for (index, count) in additions.into_iter().enumerate().filter(|(_, count)| *count > 0) {
      let c = AsciiChar::from_byte_unchecked(index as u8);
      *self.counts.entry(c).or_default() += count;
    }
  }

  pub fn add_seq_str(&mut self, sequence: &str) -> Result<(), Report> {
    self.add_seq(Self::str_to_chars(sequence)?);
    Ok(())
  }

  fn str_to_chars(s: &str) -> Result<Vec<AsciiChar>, Report> {
    s.bytes().map(AsciiChar::try_new).collect::<Result<Vec<_>, _>>()
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
