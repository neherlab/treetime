use crate::port::constants::GAP_CHAR;
use crate::port::mutation::{InDel, Mut};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct Composition {
  counts: BTreeMap<char, usize>,
}

impl Composition {
  /// Initialize counters with zeros, given an alphabet
  pub fn new<I>(alphabet_chars: I) -> Self
  where
    I: IntoIterator<Item = char>,
  {
    let mut counts: BTreeMap<char, usize> = alphabet_chars.into_iter().map(|c| (c, 0)).collect();
    counts.insert(GAP_CHAR, 0);
    Self { counts }
  }

  pub fn counts(&self) -> &BTreeMap<char, usize> {
    &self.counts
  }

  /// Initialize counters to the composition of a given sequence
  pub fn with_sequence<CI, SI>(sequence: SI, alphabet_chars: CI) -> Self
  where
    CI: IntoIterator<Item = char>,
    SI: IntoIterator<Item = char>,
  {
    let mut this = Self::new(alphabet_chars);
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

  // /// Reflect sequence mutation in the composition counts
  // pub fn add_mutation(&mut self, m: &Mut) {
  //   *self.counts.entry(m.reff).or_default() -= 1;
  //   *self.counts.entry(m.qry).or_default() += 1;
  // }
  //
  // /// Reflect sequence indel in the composition counts
  // pub fn add_indel(&mut self, indel: &InDel) {
  //   for nuc in &indel.seq {
  //     if indel.deletion {
  //       *self.counts.entry(*nuc).or_default() -= 1;
  //       *self.counts.entry(GAP_CHAR).or_default() += 1;
  //     } else {
  //       *self.counts.entry(*nuc).or_default() += 1;
  //       *self.counts.entry(GAP_CHAR).or_default() -= 1;
  //     }
  //   }
  // }

  /// Reflect sequence mutation in the composition counts
  pub fn add_mutation(&mut self, m: &Mut) {
    self.adjust_count(m.reff, -1);
    self.adjust_count(m.qry, 1);
  }

  /// Reflect sequence indel in the composition counts
  pub fn add_indel(&mut self, indel: &InDel) {
    let adjust_by = if indel.deletion { -1 } else { 1 };
    for nuc in &indel.seq {
      self.adjust_count(*nuc, adjust_by);
      self.adjust_count(GAP_CHAR, -adjust_by);
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
  use crate::port::mutation::Mut;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_composition_empty() {
    let comp = Composition::new("ACGT".chars());
    assert_eq!(
      comp.counts(),
      &btreemap! { GAP_CHAR => 0, 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0}
    );
  }

  #[test]
  fn test_composition_with_sequence() {
    let comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT".chars());
    assert_eq!(
      &btreemap! { GAP_CHAR => 0, 'A' => 6, 'C' => 5, 'G' => 6, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_sequence() {
    let mut comp = Composition::new("ACGT".chars());
    comp.add_sequence("AAAGCTTACGGGGTCAAGTCC".chars());
    assert_eq!(
      &btreemap! { GAP_CHAR => 0, 'A' => 6, 'C' => 5, 'G' => 6, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_sequence_with_refs() {
    let mut comp = Composition::new("ACGT".chars());
    let sequence: Vec<char> = "AAAGCTTACGGGGTCAAGTCC".chars().collect();
    comp.add_sequence(sequence.iter().copied());
    assert_eq!(
      &btreemap! { GAP_CHAR => 0, 'A' => 6, 'C' => 5, 'G' => 6, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_mutation() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT".chars());
    let mutation = Mut {
      pos: 123,
      reff: 'A',
      qry: 'G',
    };
    comp.add_mutation(&mutation);
    assert_eq!(
      &btreemap! { GAP_CHAR => 0, 'A' => 5, 'C' => 5, 'G' => 7, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_deletion() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT".chars());

    let indel = InDel {
      range: (1, 5),
      seq: vec!['A', 'A', 'G', 'C'],
      deletion: true,
    };
    comp.add_indel(&indel);
    assert_eq!(
      &btreemap! { GAP_CHAR => 4, 'A' => 4, 'C' => 4, 'G' => 5, 'T' => 4},
      comp.counts()
    );
  }

  #[test]
  fn test_composition_add_insertion() {
    let mut comp = Composition::with_sequence("AAAGCTTACGGGGTCAAGTCC".chars(), "ACGT".chars());
    let indel = InDel {
      range: (3, 6),
      seq: vec!['A', 'T', 'C'],
      deletion: false,
    };
    comp.add_indel(&indel);
    assert_eq!(
      &btreemap! { GAP_CHAR => 0, 'A' => 7, 'C' => 6, 'G' => 6, 'T' => 5},
      comp.counts()
    );
  }
}
