use crate::representation::seq_char::AsciiChar;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::representation::seq::Seq;
  use crate::seq::mutation::Sub;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;

  #[test]
  fn test_composition_empty() {
    let actual = Composition::new("ACGT-".bytes(), b'-');
    let expected = Composition::from(btreemap! { b'-' => 0, b'A' => 0, b'C' => 0, b'G' => 0, b'T' => 0}, b'-');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_with_sequence() {
    let actual = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    let expected = Composition::from(btreemap! { b'-' => 0, b'A' => 6, b'C' => 5, b'G' => 6, b'T' => 4}, b'-');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_with_sequence_and_alphabet() {
    let chars = Alphabet::new(AlphabetName::Nuc, false).unwrap().chars().collect_vec();
    let actual = Composition::with_sequence("ACATCGCCNNA--GAC".bytes(), chars, b'-');
    let expected = Composition::from(
      btreemap! {
        b'-' => 2,
        b'A' => 4,
        b'B' => 0,
        b'C' => 5,
        b'D' => 0,
        b'G' => 2,
        b'H' => 0,
        b'K' => 0,
        b'M' => 0,
        b'N' => 2,
        b'R' => 0,
        b'S' => 0,
        b'T' => 1,
        b'V' => 0,
        b'W' => 0,
        b'Y' => 0,
      },
      b'-',
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_sequence() {
    let mut actual = Composition::new("ACGT-".bytes(), b'-');
    actual.add_sequence("AAAGCTTACGGGGTCAAGTCC".bytes());
    let expected = Composition::from(btreemap! { b'-' => 0, b'A' => 6, b'C' => 5, b'G' => 6, b'T' => 4}, b'-');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_sequence_with_refs() {
    let mut actual = Composition::new("ACGT-".bytes(), b'-');
    let sequence: Seq = "AAAGCTTACGGGGTCAAGTCC".into();
    actual.add_sequence(sequence);
    let expected = Composition::from(btreemap! { b'-' => 0, b'A' => 6, b'C' => 5, b'G' => 6, b'T' => 4}, b'-');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_mutation() {
    let mut actual = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    let mutation = Sub::from_str("A123G").unwrap();
    actual.add_sub(&mutation);
    let expected = Composition::from(btreemap! { b'-' => 0, b'A' => 5, b'C' => 5, b'G' => 7, b'T' => 4}, b'-');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_deletion() {
    let mut actual = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    let indel = InDel::del((1, 5), "AAGC");
    actual.add_indel(&indel);
    let expected = Composition::from(btreemap! { b'-' => 4, b'A' => 4, b'C' => 4, b'G' => 5, b'T' => 4}, b'-');
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_insertion() {
    let mut actual = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    let indel = InDel::ins((3, 6), "ATC");
    actual.add_indel(&indel);
    let expected = Composition::from(btreemap! { b'-' => 0, b'A' => 7, b'C' => 6, b'G' => 6, b'T' => 5}, b'-');
    assert_eq!(expected, actual);
  }
}
