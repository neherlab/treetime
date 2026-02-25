//! Proptest generators for sequence alignments.

use proptest::prelude::*;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::seq::Seq;

/// Generate a single nucleotide character.
///
/// 90% standard (ACGT), 10% ambiguous/gaps (N, -, R, Y).
pub fn arb_nucleotide() -> impl Strategy<Value = char> {
  prop_oneof![
    9 => prop::sample::select(vec!['A', 'C', 'G', 'T']),
    1 => prop::sample::select(vec!['N', '-', 'R', 'Y']),
  ]
}

/// Generate a nucleotide sequence of given length.
pub fn arb_sequence(len: usize) -> impl Strategy<Value = String> {
  prop::collection::vec(arb_nucleotide(), len).prop_map(|chars| chars.into_iter().collect())
}

/// Generate a gap-free nucleotide (ACGT only).
pub fn arb_nucleotide_no_gaps() -> impl Strategy<Value = char> {
  prop::sample::select(vec!['A', 'C', 'G', 'T'])
}

/// Generate a gap-free sequence (no ambiguous characters).
pub fn arb_sequence_no_gaps(len: usize) -> impl Strategy<Value = String> {
  prop::collection::vec(arb_nucleotide_no_gaps(), len).prop_map(|chars| chars.into_iter().collect())
}

/// Generate an alignment (FASTA records) for given taxa names.
///
/// All sequences have the same length.
pub fn arb_alignment(taxa: Vec<String>, seq_len: usize) -> impl Strategy<Value = Vec<FastaRecord>> {
  let n = taxa.len();
  prop::collection::vec(arb_sequence(seq_len), n).prop_map(move |sequences| {
    taxa
      .iter()
      .zip(sequences)
      .map(|(name, seq_str)| FastaRecord {
        seq_name: name.clone(),
        seq: Seq::try_from_str(&seq_str).expect("Generated sequence should be valid ASCII"),
        ..FastaRecord::default()
      })
      .collect()
  })
}

/// Generate a gap-free alignment for given taxa names.
pub fn arb_alignment_no_gaps(taxa: Vec<String>, seq_len: usize) -> impl Strategy<Value = Vec<FastaRecord>> {
  let n = taxa.len();
  prop::collection::vec(arb_sequence_no_gaps(seq_len), n).prop_map(move |sequences| {
    taxa
      .iter()
      .zip(sequences)
      .map(|(name, seq_str)| FastaRecord {
        seq_name: name.clone(),
        seq: Seq::try_from_str(&seq_str).expect("Generated sequence should be valid ASCII"),
        ..FastaRecord::default()
      })
      .collect()
  })
}

#[cfg(test)]
mod tests {
  use super::*;
  use proptest::proptest;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn test_prop_alignment_arb_sequence_length(seq in arb_sequence(10)) {
      prop_assert_eq!(seq.len(), 10);
    }

    #[test]
    fn test_prop_alignment_arb_sequence_valid_chars(seq in arb_sequence(20)) {
      for c in seq.chars() {
        prop_assert!(
          ['A', 'C', 'G', 'T', 'N', '-', 'R', 'Y'].contains(&c),
          "Invalid character: {c}"
        );
      }
    }

    #[test]
    fn test_prop_alignment_arb_alignment_structure(
      aln in arb_alignment(vec!["A".to_owned(), "B".to_owned(), "C".to_owned()], 15)
    ) {
      prop_assert_eq!(aln.len(), 3);
      prop_assert_eq!(&aln[0].seq_name, "A");
      prop_assert_eq!(&aln[1].seq_name, "B");
      prop_assert_eq!(&aln[2].seq_name, "C");
      for record in &aln {
        prop_assert_eq!(record.seq.len(), 15);
      }
    }

    #[test]
    fn test_prop_alignment_arb_alignment_no_gaps_chars(
      aln in arb_alignment_no_gaps(vec!["X".to_owned(), "Y".to_owned()], 10)
    ) {
      for record in &aln {
        for c in &record.seq {
          let ch = char::from(*c);
          prop_assert!(
            ['A', 'C', 'G', 'T'].contains(&ch),
            "Gap-free alignment contains invalid char: {ch}"
          );
        }
      }
    }
  }
}
