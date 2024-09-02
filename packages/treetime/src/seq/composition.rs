use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Composition {
  counts: BTreeMap<char, usize>,
  gap: char,
}

impl Composition {
  /// Initialize counters with zeros, given an alphabet
  pub fn new<I>(alphabet_chars: I, gap: char) -> Self
  where
    I: IntoIterator<Item = char>,
  {
    let counts: BTreeMap<char, usize> = alphabet_chars.into_iter().map(|c| (c, 0)).collect();
    Self { counts, gap }
  }

  pub fn get(&self, c: char) -> Option<usize> {
    self.counts.get(&c).copied()
  }

  pub fn counts(&self) -> &BTreeMap<char, usize> {
    &self.counts
  }

  /// Initialize counters to the composition of a given sequence
  pub fn with_sequence<CI, SI>(sequence: SI, alphabet_chars: CI, gap: char) -> Self
  where
    CI: IntoIterator<Item = char>,
    SI: IntoIterator<Item = char>,
  {
    let mut this = Self::new(alphabet_chars, gap);
    this.add_sequence(sequence);
    this
  }

  /// Add composition of a given sequence to the counts
  pub fn add_sequence<SI>(&mut self, sequence: SI)
  where
    SI: IntoIterator<Item = char>,
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

  pub fn adjust_count(&mut self, nuc: char, change: isize) {
    let count = self.counts.entry(nuc).or_default();
    *count = count.saturating_add_signed(change);
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::mutation::Sub;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_composition_empty() {
    let comp = Composition::new("ACGT-".chars(), '-');
    assert_eq!(
      comp.counts(),
      &btreemap! { '-' => 0, 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0}
    );
  }

  #[test]
  fn test_composition_with_sequence() {
    let comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT-".chars(), '-');
    assert_eq!(
      &btreemap! { '-' => 0, 'A' => 6, 'C' => 5, 'G' => 6, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_with_sequence_and_alphabet() {
    let chars = Alphabet::new(AlphabetName::Nuc, false).unwrap().chars().collect_vec();
    let comp = Composition::with_sequence("ACATCGCCNNA--GAC".chars(), chars, '-');
    assert_eq!(
      &btreemap! {
          '-' => 2,
          'A' => 4,
          'B' => 0,
          'C' => 5,
          'D' => 0,
          'G' => 2,
          'H' => 0,
          'K' => 0,
          'M' => 0,
          'N' => 2,
          'R' => 0,
          'S' => 0,
          'T' => 1,
          'V' => 0,
          'W' => 0,
          'Y' => 0,
      },
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_sequence() {
    let mut comp = Composition::new("ACGT-".chars(), '-');
    comp.add_sequence("AAAGCTTACGGGGTCAAGTCC".chars());
    assert_eq!(
      &btreemap! { '-' => 0, 'A' => 6, 'C' => 5, 'G' => 6, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_sequence_with_refs() {
    let mut comp = Composition::new("ACGT-".chars(), '-');
    let sequence: Vec<char> = "AAAGCTTACGGGGTCAAGTCC".chars().collect();
    comp.add_sequence(sequence.iter().copied());
    assert_eq!(
      &btreemap! { '-' => 0, 'A' => 6, 'C' => 5, 'G' => 6, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_mutation() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT-".chars(), '-');
    let mutation = Sub {
      pos: 123,
      reff: 'A',
      qry: 'G',
    };
    comp.add_sub(&mutation);
    assert_eq!(
      &btreemap! { '-' => 0, 'A' => 5, 'C' => 5, 'G' => 7, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_deletion() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT-".chars(), '-');

    let indel = InDel {
      range: (1, 5),
      seq: vec!['A', 'A', 'G', 'C'],
      deletion: true,
    };
    comp.add_indel(&indel);
    assert_eq!(
      &btreemap! { '-' => 4, 'A' => 4, 'C' => 4, 'G' => 5, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_insertion() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT-".chars(), '-');
    let indel = InDel {
      range: (3, 6),
      seq: vec!['A', 'T', 'C'],
      deletion: false,
    };
    comp.add_indel(&indel);
    assert_eq!(
      &btreemap! { '-' => 0, 'A' => 7, 'C' => 6, 'G' => 6, 'T' => 5},
      comp.counts()
    );
  }
}
