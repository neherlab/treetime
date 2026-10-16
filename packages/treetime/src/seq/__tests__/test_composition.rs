#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::composition::Composition;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{deletion, insertion};
  use eyre::Report;
  use itertools::Itertools;
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
    let composition = Composition::new(chars("ACGT-"), c(b'-'));
    let actual = chars("ACGT-N").map(|ch| composition.get(ch)).collect_vec();
    assert_eq!(vec![Some(0), Some(0), Some(0), Some(0), Some(0), None], actual);
  }

  #[test]
  fn test_composition_with_seq() -> Result<(), Report> {
    let actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, chars("ACGT-"), c(b'-'));
    let expected = helpers::composition(chars("ACGT-"), &[(b'A', 6), (b'C', 5), (b'G', 6), (b'T', 4)]);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_with_seq_and_alphabet() -> Result<(), Report> {
    let alpha_chars = Alphabet::new(AlphabetName::Nuc).unwrap().chars().collect_vec();
    let actual = Composition::with_seq(Seq::try_from_str("ACATCGCCNNA--GAC")?, alpha_chars.clone(), c(b'-'));
    let expected = helpers::composition(
      alpha_chars,
      &[(b'-', 2), (b'A', 4), (b'C', 5), (b'G', 2), (b'N', 2), (b'T', 1)],
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_add_seq_with_refs() -> Result<(), Report> {
    let mut actual = Composition::new(chars("ACGT-"), c(b'-'));
    let sequence = Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?;
    actual.add_seq(&sequence);
    let expected = helpers::composition(chars("ACGT-"), &[(b'A', 6), (b'C', 5), (b'G', 6), (b'T', 4)]);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_add_mutation() -> Result<(), Report> {
    let mut actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, chars("ACGT-"), c(b'-'));
    let mutation = Sub::from_str("A123G").unwrap();
    actual.add_sub(&mutation);
    let expected = helpers::composition(chars("ACGT-"), &[(b'A', 5), (b'C', 5), (b'G', 7), (b'T', 4)]);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_add_deletion() -> Result<(), Report> {
    let mut actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, chars("ACGT-"), c(b'-'));
    let indel = deletion((1, 5), Seq::try_from_str("AAGC")?);
    actual.add_indel(&indel);
    let expected = helpers::composition(chars("ACGT-"), &[(b'-', 4), (b'A', 4), (b'C', 4), (b'G', 5), (b'T', 4)]);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_composition_add_insertion() -> Result<(), Report> {
    let mut actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, chars("ACGT-"), c(b'-'));
    let indel = insertion((3, 6), Seq::try_from_str("ATC")?);
    actual.add_indel(&indel);
    let expected = helpers::composition(chars("ACGT-"), &[(b'A', 7), (b'C', 6), (b'G', 6), (b'T', 5)]);
    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use crate::seq::composition::Composition;
    use treetime_primitives::AsciiChar;

    pub(super) fn composition(
      alphabet_chars: impl IntoIterator<Item = AsciiChar>,
      counts: &[(u8, isize)],
    ) -> Composition {
      let mut composition = Composition::new(alphabet_chars, AsciiChar::from_byte_unchecked(b'-'));
      for &(c, count) in counts {
        composition.adjust_count(AsciiChar::from_byte_unchecked(c), count);
      }
      composition
    }
  }
}
