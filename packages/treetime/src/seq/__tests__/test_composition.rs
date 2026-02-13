#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use treetime_primitives::Seq;
  use crate::seq::composition::Composition;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use itertools::Itertools;
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
