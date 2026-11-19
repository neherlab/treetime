#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::mask::{create_mask, mask_to_string};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::Seq;

  #[rustfmt::skip]
  #[rstest]
  #[case::all_unknown_columns_masked(AlphabetName::Nuc, &["NNNN", "NNNN"], 4, vec![true, true, true, true])]
  #[case::one_observed_state_clears_column(AlphabetName::Nuc, &["NNNN", "NNNN", "ANGN"], 4, vec![false, true, false, true])]
  #[case::gaps_are_observed_states(AlphabetName::Nuc, &["N-N", "NN-", "N-N"], 3, vec![true, false, false])]
  #[case::partial_ambiguity_codes_are_observed_states(AlphabetName::Nuc, &["RNY", "NNN"], 3, vec![false, true, false])]
  #[case::amino_acid_unknown_is_x(AlphabetName::Aa, &["XM-X", "XXX-"], 4, vec![true, false, false, false])]
  #[case::empty_alignment_masks_every_column(AlphabetName::Nuc, &[], 5, vec![true, true, true, true, true])]
  #[trace]
  fn test_mask_marks_columns_without_observed_state(
    #[case] alphabet: AlphabetName,
    #[case] sequences: &[&str],
    #[case] alignment_length: usize,
    #[case] expected: Vec<bool>,
  ) {
    let alphabet = Alphabet::new(alphabet).unwrap();
    let seqs: Vec<Seq> = sequences.iter().map(|seq| Seq::try_from_str(seq).unwrap()).collect();

    let mask = create_mask(&seqs, alignment_length, &alphabet);
    assert_eq!(expected, mask);
  }

  #[test]
  fn test_mask_to_string_writes_masked_as_one() {
    let mask = vec![true, false, false, true, false];
    assert_eq!("10010", mask_to_string(&mask));
  }
}
