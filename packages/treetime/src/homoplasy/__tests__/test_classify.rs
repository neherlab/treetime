#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::homoplasy::classify::{MutationClass, classify_mutation};
  use crate::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
  use crate::test_utils::{deletion, insertion};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::str::FromStr;
  use treetime_primitives::Seq;

  #[rustfmt::skip]
  #[rstest]
  #[case::determined(           "A10G", MutationClass::Substitution)]
  #[case::to_partial_ambiguity( "G10R", MutationClass::Ambiguous)]
  #[case::to_two_state_code(    "A10Y", MutationClass::Ambiguous)]
  #[case::to_unknown(           "A10N", MutationClass::Ambiguous)]
  #[case::from_unknown(         "N10G", MutationClass::Ambiguous)]
  #[trace]
  fn test_classify_substitution(#[case] text: &str, #[case] expected: MutationClass) -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let event = MutationEvent::Substitution(Sub::from_str(text)?);
    assert_eq!(expected, classify_mutation(&event, &alphabet));
    Ok(())
  }

  #[test]
  fn test_classify_stop_codon_substitution_in_amino_acids() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Aa)?;
    let event = MutationEvent::Substitution(Sub::from_str("W5*")?);
    assert_eq!(MutationClass::Substitution, classify_mutation(&event, &alphabet));
    Ok(())
  }

  #[test]
  fn test_classify_insertion_and_deletion() -> Result<(), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let inserted = Mutation::indel(MutationTrack::Nucleotide, &insertion((3, 5), Seq::try_from_slice(b"AC")?))?;
    let deleted = Mutation::indel(MutationTrack::Nucleotide, &deletion((3, 5), Seq::try_from_slice(b"AC")?))?;
    assert_eq!(
      (MutationClass::Indel, MutationClass::Indel),
      (
        classify_mutation(&inserted.event, &alphabet),
        classify_mutation(&deleted.event, &alphabet)
      )
    );
    Ok(())
  }
}
