#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::mask::{create_mask, mask_to_string};
  use treetime_primitives::{AlignmentRecord, Seq};

  fn nuc_alphabet() -> Alphabet {
    Alphabet::new(AlphabetName::Nuc).unwrap()
  }

  fn record(name: &str, sequence: &str) -> AlignmentRecord {
    AlignmentRecord {
      name: name.to_owned(),
      seq: Seq::try_from_str(sequence).unwrap(),
    }
  }

  #[test]
  fn test_mask_no_ambiguity() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "ACGT"), record("B", "ACGT"), record("C", "TGCA")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 4, &alphabet);
    assert_eq!(vec![false, false, false, false], mask);
    assert_eq!("0000", mask_to_string(&mask));
  }

  #[test]
  fn test_mask_all_ambiguous() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "NNNN"), record("B", "NNNN")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 4, &alphabet);
    assert_eq!(vec![true, true, true, true], mask);
    assert_eq!("1111", mask_to_string(&mask));
  }

  #[test]
  fn test_mask_mixed() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "ANNG"), record("B", "NNNA"), record("C", "NNN-")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 4, &alphabet);
    assert_eq!(vec![false, true, true, false], mask);
    assert_eq!("0110", mask_to_string(&mask));
  }

  #[test]
  fn test_mask_single_informative_tip_clears_position() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "NNNN"), record("B", "NNNN"), record("C", "ANGN")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 4, &alphabet);
    assert_eq!(vec![false, true, false, true], mask);
  }

  #[test]
  fn test_mask_gaps_are_observed_states() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "----"), record("B", "----")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 4, &alphabet);
    assert_eq!(vec![false, false, false, false], mask);
  }

  #[test]
  fn test_mask_column_of_unknowns_and_gaps_is_not_masked() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "N-N"), record("B", "NN-"), record("C", "N-N")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 3, &alphabet);
    assert_eq!(vec![true, false, false], mask);
  }

  #[test]
  fn test_mask_partial_ambiguity_codes_are_not_masked() {
    let alphabet = nuc_alphabet();
    let aln = [record("A", "RNY"), record("B", "NNN")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 3, &alphabet);
    assert_eq!(vec![false, true, false], mask);
  }

  #[test]
  fn test_mask_amino_acids_uses_x_as_unknown() {
    let alphabet = Alphabet::new(AlphabetName::Aa).unwrap();
    let aln = [record("A", "XM-X"), record("B", "XXX-")];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 4, &alphabet);
    assert_eq!(vec![true, false, false, false], mask);
  }

  #[test]
  fn test_mask_empty_alignment() {
    let alphabet = nuc_alphabet();
    let aln: Vec<AlignmentRecord> = vec![];

    let mask = create_mask(aln.iter().map(|record| &record.seq), 5, &alphabet);
    assert_eq!(vec![true, true, true, true, true], mask);
  }

  #[test]
  fn test_mask_to_string_roundtrip() {
    let mask = vec![true, false, false, true, false];
    assert_eq!("10010", mask_to_string(&mask));
  }
}
