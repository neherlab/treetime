#[cfg(test)]
mod tests {
  use self::helpers::{c, nonzero_counts};
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::seq::composition::Composition;
  use crate::seq::mutation::Sub;
  use crate::test_utils::{deletion, insertion};
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::str::FromStr;
  use treetime_primitives::{AlphabetLike, Seq};

  #[test]
  fn test_composition_empty_counts_zero_for_every_character_of_the_alphabet() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let composition = Composition::new(&alphabet);

    let counts = alphabet.chars().map(|ch| composition.get(ch)).collect_vec();
    assert_eq!(vec![Some(0); alphabet.chars().count()], counts);
    Ok(())
  }

  #[test]
  fn test_composition_has_no_count_for_a_character_outside_the_alphabet() -> Result<(), Report> {
    let composition = Composition::new(&Alphabet::new(AlphabetName::Nuc)?);

    assert_eq!(None, composition.get(c(b'*')));
    Ok(())
  }

  #[test]
  fn test_composition_with_seq_counts_every_character() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;

    let actual = Composition::with_seq(Seq::try_from_str("ACATCGCCNNA--GAC")?, &alphabet);

    let expected = vec![('-', 2), ('A', 4), ('C', 5), ('G', 2), ('N', 2), ('T', 1)];
    assert_eq!(expected, nonzero_counts(&actual, &alphabet));
    Ok(())
  }

  #[test]
  fn test_composition_add_seq_adds_to_the_counts() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let mut actual = Composition::with_seq(Seq::try_from_str("ACGT")?, &alphabet);

    actual.add_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?);

    let expected = vec![('A', 7), ('C', 6), ('G', 7), ('T', 5)];
    assert_eq!(expected, nonzero_counts(&actual, &alphabet));
    Ok(())
  }

  #[test]
  fn test_composition_add_sub_moves_one_count() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let mut actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, &alphabet);

    actual.add_sub(&Sub::from_str("A123G")?);

    let expected = vec![('A', 5), ('C', 5), ('G', 7), ('T', 4)];
    assert_eq!(expected, nonzero_counts(&actual, &alphabet));
    Ok(())
  }

  #[test]
  fn test_composition_add_deletion_moves_counts_to_the_gap() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let mut actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, &alphabet);

    actual.add_indel(&deletion((1, 5), Seq::try_from_str("AAGC")?));

    let expected = vec![('-', 4), ('A', 4), ('C', 4), ('G', 5), ('T', 4)];
    assert_eq!(expected, nonzero_counts(&actual, &alphabet));
    Ok(())
  }

  #[test]
  fn test_composition_add_insertion_keeps_the_gap_count_at_zero() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let mut actual = Composition::with_seq(Seq::try_from_str("AAAGCTTACGGGGTCAAGTCC")?, &alphabet);

    actual.add_indel(&insertion((3, 6), Seq::try_from_str("ATC")?));

    let expected = vec![('A', 7), ('C', 6), ('G', 6), ('T', 5)];
    assert_eq!(expected, nonzero_counts(&actual, &alphabet));
    Ok(())
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::seq::composition::Composition;
    use treetime_primitives::{AlphabetLike, AsciiChar};

    pub(super) fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub(super) fn nonzero_counts(composition: &Composition, alphabet: &Alphabet) -> Vec<(char, usize)> {
      alphabet
        .chars()
        .filter_map(|ch| {
          composition
            .get(ch)
            .filter(|&count| count > 0)
            .map(|count| (char::from(ch), count))
        })
        .collect()
    }
  }
}
