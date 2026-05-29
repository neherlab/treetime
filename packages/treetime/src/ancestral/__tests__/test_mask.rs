#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::mask::{create_mask, mask_to_string};
  use treetime_io::fasta::FastaRecord;
  use treetime_primitives::Seq;

  fn nuc_alphabet() -> Alphabet {
    Alphabet::new(AlphabetName::Nuc).unwrap()
  }

  fn record(name: &str, sequence: &str) -> FastaRecord {
    FastaRecord {
      seq_name: name.to_owned(),
      seq: Seq::try_from_str(sequence).unwrap(),
      ..Default::default()
    }
  }

  #[test]
  fn test_mask_no_ambiguity() {
    let alphabet = nuc_alphabet();
    let aln = vec![record("A", "ACGT"), record("B", "ACGT"), record("C", "TGCA")];

    let mask = create_mask(&aln, 4, &alphabet);
    assert_eq!(vec![false, false, false, false], mask);
    assert_eq!("0000", mask_to_string(&mask));
  }

  #[test]
  fn test_mask_all_ambiguous() {
    let alphabet = nuc_alphabet();
    let aln = vec![record("A", "NNNN"), record("B", "NNNN")];

    let mask = create_mask(&aln, 4, &alphabet);
    assert_eq!(vec![true, true, true, true], mask);
    assert_eq!("1111", mask_to_string(&mask));
  }

  #[test]
  fn test_mask_mixed() {
    let alphabet = nuc_alphabet();
    let aln = vec![record("A", "ANNG"), record("B", "NNNA"), record("C", "NNN-")];

    let mask = create_mask(&aln, 4, &alphabet);
    assert_eq!(vec![false, true, true, false], mask);
    assert_eq!("0110", mask_to_string(&mask));
  }

  #[test]
  fn test_mask_single_informative_tip_clears_position() {
    let alphabet = nuc_alphabet();
    let aln = vec![record("A", "NNNN"), record("B", "NNNN"), record("C", "ANGN")];

    let mask = create_mask(&aln, 4, &alphabet);
    assert_eq!(vec![false, true, false, true], mask);
  }

  #[test]
  fn test_mask_gaps_treated_as_ambiguous() {
    let alphabet = nuc_alphabet();
    let aln = vec![record("A", "----"), record("B", "----")];

    let mask = create_mask(&aln, 4, &alphabet);
    assert_eq!(vec![true, true, true, true], mask);
  }

  #[test]
  fn test_mask_empty_alignment() {
    let alphabet = nuc_alphabet();
    let aln: Vec<FastaRecord> = vec![];

    let mask = create_mask(&aln, 5, &alphabet);
    assert_eq!(vec![true, true, true, true, true], mask);
  }

  #[test]
  fn test_mask_to_string_roundtrip() {
    let mask = vec![true, false, false, true, false];
    assert_eq!("10010", mask_to_string(&mask));
  }
}
