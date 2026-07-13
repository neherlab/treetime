#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::__tests__::prop_generators::alignment::arb_alignment_no_gaps;
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::partition::fitch::PartitionFitch;
  use crate::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use proptest::prelude::*;
  use std::collections::BTreeMap;
  use std::iter;
  use treetime_graph::node::Named;
  use treetime_io::fasta::FastaRecord;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AsciiChar, Seq};

  const BALANCED: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root:0.0;";
  const ROOTED_AT_A: &str = "(A:0.05,(B:0.1,(C:0.1,D:0.1)CD:0.2)AB:0.05)root:0.0;";

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    /// Fitch 1971: on a binary tree, the backward/forward assignment attains
    /// the minimum number of unordered character-state changes.
    #[test]
    fn test_prop_fitch_informative_score_matches_exhaustive_assignment(
      alignment in arb_four_taxon_alignment(),
    ) {
      let actual = run_fitch(BALANCED, &alignment).unwrap();
      let expected = exhaustive_balanced_score(&alignment);

      prop_assert_eq!(expected, actual.score);
    }

    /// The minimum unordered-character score belongs to the unrooted binary
    /// topology, so inserting the degree-two root on another edge preserves it.
    #[test]
    fn test_prop_fitch_informative_score_preserved_by_edge_reroot(
      alignment in arb_four_taxon_alignment(),
    ) {
      let balanced = run_fitch(BALANCED, &alignment).unwrap();
      let rerooted = run_fitch(ROOTED_AT_A, &alignment).unwrap();

      prop_assert_eq!(balanced.score, rerooted.score);
    }

    #[test]
    fn test_prop_fitch_informative_column_reversal_reverses_reconstruction(
      alignment in arb_four_taxon_alignment(),
    ) {
      let original = run_fitch(BALANCED, &alignment).unwrap();
      let reversed_alignment = alignment
        .iter()
        .cloned()
        .map(|mut record| {
          record.seq = record.seq.iter().copied().rev().collect();
          record
        })
        .collect::<Vec<_>>();
      let reversed = run_fitch(BALANCED, &reversed_alignment).unwrap();
      let expected_sequences: BTreeMap<String, Seq> = original
        .sequences
        .into_iter()
        .map(|(name, sequence)| (name, sequence.iter().copied().rev().collect::<Seq>()))
        .collect();

      prop_assert_eq!(original.score, reversed.score);
      prop_assert_eq!(expected_sequences, reversed.sequences);
    }

    /// A column containing the same canonical state at every leaf has zero
    /// parsimony cost and cannot alter assignments at independent columns.
    #[test]
    fn test_prop_fitch_informative_invariant_column_preserves_reconstruction(
      alignment in arb_four_taxon_alignment(),
    ) {
      let original = run_fitch(BALANCED, &alignment).unwrap();
      let with_invariant = alignment
        .iter()
        .cloned()
        .map(|mut record| {
          record.seq = iter::once(AsciiChar::from_byte_unchecked(b'A'))
            .chain(record.seq.iter().copied())
            .collect();
          record
        })
        .collect::<Vec<_>>();
      let actual = run_fitch(BALANCED, &with_invariant).unwrap();

      prop_assert_eq!(original.score, actual.score);
      for (name, expected) in original.sequences {
        let sequence = &actual.sequences[&name];
        prop_assert_eq!(Some(&AsciiChar::from_byte_unchecked(b'A')), sequence.first());
        prop_assert_eq!(expected.as_slice(), &sequence.as_slice()[1..]);
      }
    }
  }

  fn arb_four_taxon_alignment() -> impl Strategy<Value = Vec<FastaRecord>> {
    (1_usize..=8).prop_flat_map(|length| {
      arb_alignment_no_gaps(
        vec!["A".to_owned(), "B".to_owned(), "C".to_owned(), "D".to_owned()],
        length,
      )
    })
  }

  fn run_fitch(newick: &str, alignment: &[FastaRecord]) -> Result<FitchResult, Report> {
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let partition = create_fitch_partition(&graph, 0, Alphabet::default(), alignment)?;
    let score = partition.edges.values().map(|edge| edge.fitch_subs().len()).sum();
    let sequences = collect_sequences(&graph, &partition);
    Ok(FitchResult { score, sequences })
  }

  struct FitchResult {
    score: usize,
    sequences: BTreeMap<String, Seq>,
  }

  fn collect_sequences(graph: &GraphAncestral, partition: &PartitionFitch) -> BTreeMap<String, Seq> {
    graph
      .get_nodes()
      .iter()
      .map(|node| {
        let node = node.read_arc();
        let name = node.payload().read_arc().name().unwrap().as_ref().to_owned();
        (name, partition.nodes[&node.key()].seq.sequence.clone())
      })
      .collect()
  }

  fn exhaustive_balanced_score(alignment: &[FastaRecord]) -> usize {
    let leaves = alignment
      .iter()
      .map(|record| (record.seq_name.as_str(), &record.seq))
      .collect::<BTreeMap<_, _>>();
    (0..alignment[0].seq.len())
      .map(|pos| exhaustive_balanced_column([leaves["A"][pos], leaves["B"][pos], leaves["C"][pos], leaves["D"][pos]]))
      .sum()
  }

  fn exhaustive_balanced_column(leaves: [AsciiChar; 4]) -> usize {
    let canonical = [b'A', b'C', b'G', b'T'].map(AsciiChar::from_byte_unchecked);
    canonical
      .iter()
      .flat_map(|&root| canonical.iter().map(move |&left| (root, left)))
      .flat_map(|(root, left)| canonical.iter().map(move |&right| (root, left, right)))
      .map(|(root, left, right)| {
        usize::from(root != left)
          + usize::from(root != right)
          + usize::from(left != leaves[0])
          + usize::from(left != leaves[1])
          + usize::from(right != leaves[2])
          + usize::from(right != leaves[3])
      })
      .min()
      .unwrap()
  }
}
