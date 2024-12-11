use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Composition {
  counts: BTreeMap<u8, usize>,
  gap: u8,
}

impl Composition {
  /// Initialize counters with zeros, given an alphabet
  pub fn new<I, C>(alphabet_chars: I, gap: u8) -> Self
  where
    I: IntoIterator<Item = C>,
    C: Into<u8>,
  {
    let counts: BTreeMap<u8, usize> = alphabet_chars.into_iter().map(|c| (c.into(), 0)).collect();
    Self { counts, gap }
  }

  pub fn get(&self, c: u8) -> Option<usize> {
    self.counts.get(&c).copied()
  }

  pub fn counts(&self) -> &BTreeMap<u8, usize> {
    &self.counts
  }

  /// Initialize counters to the composition of a given sequence
  pub fn with_sequence<CI, SI>(sequence: SI, alphabet_chars: CI, gap: u8) -> Self
  where
    CI: IntoIterator<Item = u8>,
    SI: IntoIterator<Item = u8>,
  {
    let mut this = Self::new(alphabet_chars, gap);
    this.add_sequence(sequence);
    this
  }

  /// Add composition of a given sequence to the counts
  pub fn add_sequence<SI>(&mut self, sequence: SI)
  where
    SI: IntoIterator<Item = u8>,
  {
    let counts = sequence.into_iter().counts();
    self.counts.extend(counts);
  }

  /// Reflect sequence mutation in the composition counts
  pub fn add_sub(&mut self, sub: &Sub) {
    self.adjust_count(sub.reff, -1);
    self.adjust_count(sub.qry, 1);
  }

  /// Reflect sequence indel in the composition counts
  pub fn add_indel(&mut self, indel: &InDel) {
    let adjust_by = if indel.deletion { -1 } else { 1 };
    for nuc in &indel.seq {
      self.adjust_count(*nuc, adjust_by);
      self.adjust_count(self.gap, -adjust_by);
    }
  }

  pub fn adjust_count(&mut self, nuc: u8, change: isize) {
    let count = self.counts.entry(nuc).or_default();
    *count = count.saturating_add_signed(change);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::representation::seq::Seq;
  use crate::seq;
  use crate::seq::mutation::Sub;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_composition_empty() {
    let comp = Composition::new("ACGT-".bytes(), b'-');
    assert_eq!(
      comp.counts(),
      &btreemap! { b'-' => 0, b'A' => 0, b'C' => 0, b'G' => 0, b'T' => 0}
    );
  }

  #[test]
  fn test_composition_with_sequence() {
    let comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    assert_eq!(
      &btreemap! { b'-' => 0, b'A' => 6, b'C' => 5, b'G' => 6, b'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_with_sequence_and_alphabet() {
    let chars = Alphabet::new(AlphabetName::Nuc, false).unwrap().chars().collect_vec();
    let comp = Composition::with_sequence("ACATCGCCNNA--GAC".bytes(), chars, b'-');
    assert_eq!(
      &btreemap! {
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
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_sequence() {
    let mut comp = Composition::new("ACGT-".bytes(), b'-');
    comp.add_sequence("AAAGCTTACGGGGTCAAGTCC".bytes());
    assert_eq!(
      &btreemap! { b'-' => 0, b'A' => 6, b'C' => 5, b'G' => 6, b'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_sequence_with_refs() {
    let mut comp = Composition::new("ACGT-".bytes(), b'-');
    let sequence: Seq = "AAAGCTTACGGGGTCAAGTCC".into();
    comp.add_sequence(sequence);
    assert_eq!(
      &btreemap! { b'-' => 0, b'A' => 6, b'C' => 5, b'G' => 6, b'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_mutation() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    let mutation = Sub {
      pos: 123,
      reff: b'A',
      qry: b'G',
    };
    comp.add_sub(&mutation);
    assert_eq!(
      &btreemap! { b'-' => 0, b'A' => 5, b'C' => 5, b'G' => 7, b'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_deletion() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');

    let indel = InDel {
      range: (1, 5),
      seq: seq!['A', 'A', 'G', 'C'],
      deletion: true,
    };
    comp.add_indel(&indel);
    assert_eq!(
      &btreemap! { b'-' => 4, b'A' => 4, b'C' => 4, b'G' => 5, b'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_insertion() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".bytes(), "ACGT-".bytes(), b'-');
    let indel = InDel {
      range: (3, 6),
      seq: seq!['A', 'T', 'C'],
      deletion: false,
    };
    comp.add_indel(&indel);
    assert_eq!(
      &btreemap! { b'-' => 0, b'A' => 7, b'C' => 6, b'G' => 6, b'T' => 5},
      comp.counts()
    );
  }
}
