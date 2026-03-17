#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::composition::Composition;
  use crate::seq::indel::InDel;
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;
  use treetime_primitives::{AlphabetLike, AsciiChar, Seq};

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn chars(s: &str) -> impl Iterator<Item = AsciiChar> + '_ {
    s.bytes().map(AsciiChar::from_byte_unchecked)
  }

  #[test]
  fn test_composition_empty() {
    let actual = Composition::new(chars("ACGT-"), c(b'-'));
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 0, c(b'A') => 0, c(b'C') => 0, c(b'G') => 0, c(b'T') => 0},
      c(b'-'),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_with_sequence() {
    let actual = Composition::with_sequence(chars("AAAGCTTACGGGGTCAAGTCC"), chars("ACGT-"), c(b'-'));
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 0, c(b'A') => 6, c(b'C') => 5, c(b'G') => 6, c(b'T') => 4},
      c(b'-'),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_with_sequence_and_alphabet() {
    let alpha_chars = Alphabet::new(AlphabetName::Nuc).unwrap().chars().collect_vec();
    let actual = Composition::with_sequence(chars("ACATCGCCNNA--GAC"), alpha_chars, c(b'-'));
    let expected = Composition::from_counts(
      btreemap! {
        c(b'-') => 2,
        c(b'A') => 4,
        c(b'B') => 0,
        c(b'C') => 5,
        c(b'D') => 0,
        c(b'G') => 2,
        c(b'H') => 0,
        c(b'K') => 0,
        c(b'M') => 0,
        c(b'N') => 2,
        c(b'R') => 0,
        c(b'S') => 0,
        c(b'T') => 1,
        c(b'V') => 0,
        c(b'W') => 0,
        c(b'Y') => 0,
      },
      c(b'-'),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_sequence() {
    let mut actual = Composition::new(chars("ACGT-"), c(b'-'));
    actual.add_sequence(chars("AAAGCTTACGGGGTCAAGTCC"));
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 0, c(b'A') => 6, c(b'C') => 5, c(b'G') => 6, c(b'T') => 4},
      c(b'-'),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_sequence_with_refs() -> Result<(), Report> {
    let mut actual = Composition::new(chars("ACGT-"), c(b'-'));
    let sequence = Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?;
    actual.add_sequence(sequence);
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 0, c(b'A') => 6, c(b'C') => 5, c(b'G') => 6, c(b'T') => 4},
      c(b'-'),
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_add_mutation() {
    let mut actual = Composition::with_sequence(chars("AAAGCTTACGGGGTCAAGTCC"), chars("ACGT-"), c(b'-'));
    let mutation = Sub::from_str("A123G").unwrap();
    actual.add_sub(&mutation);
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 0, c(b'A') => 5, c(b'C') => 5, c(b'G') => 7, c(b'T') => 4},
      c(b'-'),
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_composition_add_deletion() -> Result<(), Report> {
    let mut actual = Composition::with_sequence(chars("AAAGCTTACGGGGTCAAGTCC"), chars("ACGT-"), c(b'-'));
    let indel = InDel::del((1, 5), Seq::try_from_str("AAGC")?);
    actual.add_indel(&indel);
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 4, c(b'A') => 4, c(b'C') => 4, c(b'G') => 5, c(b'T') => 4},
      c(b'-'),
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_add_insertion() -> Result<(), Report> {
    let mut actual = Composition::with_sequence(chars("AAAGCTTACGGGGTCAAGTCC"), chars("ACGT-"), c(b'-'));
    let indel = InDel::ins((3, 6), Seq::try_from_str("ATC")?);
    actual.add_indel(&indel);
    let expected = Composition::from_counts(
      btreemap! { c(b'-') => 0, c(b'A') => 7, c(b'C') => 6, c(b'G') => 6, c(b'T') => 5},
      c(b'-'),
    );
    assert_eq!(expected, actual);
    Ok(())
  }
}
