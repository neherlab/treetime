#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::attach::{MissingLeafLimit, complete_alignment_for_leaves, sanitize_to_alphabet};
  use crate::ancestral::mask::mask_to_string;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::pair_leaf_sequences;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node_names::NodeNames;

  use treetime_io::nwk::nwk_read;
  use treetime_primitives::{AlignmentRecord, Seq};
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_attach_synthesizes_all_unknown_for_missing_tip() {
    let (graph, names) = helpers::four_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT"), ("C", "ACGT")]);

    let by_name = helpers::complete(
      &graph,
      &names,
      sequences,
      &alphabet,
      MissingLeafLimit { ignored: false },
    )
    .unwrap();

    assert_eq!(4, by_name.len());
    assert_eq!(Seq::try_from_str("ACGT").unwrap(), by_name["A"]);
    assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["D"]);
  }

  #[test]
  fn test_attach_aborts_above_one_third_missing() {
    let (graph, names) = helpers::two_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGT")]);

    assert_error!(
      helpers::complete(
        &graph,
        &names,
        sequences,
        &alphabet,
        MissingLeafLimit { ignored: false }
      ),
      "More than one third of terminal nodes (1 of 2) cannot be assigned a sequence. \
       Are you sure the alignment belongs to the tree? \
       Pass --allow-missing-sequences to proceed, treating missing tips as fully ambiguous."
    );
  }

  #[test]
  fn test_attach_allow_missing_sequences_bypasses_threshold() {
    let (graph, names) = helpers::two_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGT")]);

    let by_name = helpers::complete(&graph, &names, sequences, &alphabet, MissingLeafLimit { ignored: true }).unwrap();

    assert_eq!(2, by_name.len());
    assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["B"]);
  }

  #[test]
  fn test_attach_exactly_one_third_missing_does_not_abort() {
    let (graph, names) = helpers::three_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT")]);

    let by_name = helpers::complete(
      &graph,
      &names,
      sequences,
      &alphabet,
      MissingLeafLimit { ignored: false },
    )
    .unwrap();

    assert_eq!(3, by_name.len());
    assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["C"]);
  }

  #[test]
  fn test_attach_keeps_extra_records_not_matching_any_leaf() {
    let (graph, names) = helpers::two_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT"), ("reference", "ACGT")]);
    let mut paired = pair_leaf_sequences(&graph, &names, sequences).sequences;

    complete_alignment_for_leaves(
      &graph,
      &mut paired.nodes,
      4,
      &alphabet,
      MissingLeafLimit { ignored: false },
      "alignment 'aln.fasta'",
      &NoopProgress,
    )
    .unwrap();

    assert_eq!(
      vec![(o!("reference"), Seq::try_from_str("ACGT").unwrap())],
      paired.unmatched
    );
  }

  #[test]
  fn test_attach_mask_ignores_records_not_matching_any_leaf() {
    let (graph, names) = helpers::three_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACNN-"), ("B", "AGNN-"), ("C", "ATNN-"), ("outlier", "AAGTG")]);
    let mut paired = pair_leaf_sequences(&graph, &names, sequences).sequences;
    let length = paired.common_length().unwrap();
    complete_alignment_for_leaves(
      &graph,
      &mut paired.nodes,
      length,
      &alphabet,
      MissingLeafLimit { ignored: false },
      "alignment 'aln.fasta'",
      &NoopProgress,
    )
    .unwrap();

    let mask = paired.mask(length, &alphabet);

    assert_eq!("00110", mask_to_string(&mask));
  }

  #[test]
  fn test_attach_mask_counts_leaves_without_a_record_as_unknown() {
    let (graph, names) = helpers::three_leaf_tree();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGN"), ("B", "NNGN"), ("outlier", "AAAA")]);
    let mut paired = pair_leaf_sequences(&graph, &names, sequences).sequences;
    let length = paired.common_length().unwrap();
    complete_alignment_for_leaves(
      &graph,
      &mut paired.nodes,
      length,
      &alphabet,
      MissingLeafLimit { ignored: true },
      "alignment 'aln.fasta'",
      &NoopProgress,
    )
    .unwrap();

    let mask = paired.mask(length, &alphabet);

    assert_eq!("0001", mask_to_string(&mask));
  }

  #[test]
  fn test_attach_gives_leaves_that_share_a_name_the_first_record_of_that_name() {
    let nwk_parsed = nwk_read(b"(A:0.1,A:0.1,B:0.1)root;".as_slice()).unwrap();
    let names = nwk_parsed.names.clone();
    let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
    let sequences = helpers::records(&[("A", "ACGT"), ("B", "GGGG"), ("A", "TTTT")]);

    let pairing = pair_leaf_sequences(&nwk_parsed.graph, &names, sequences);
    let leaves = nwk_parsed.graph.get_leaves().map(|leaf| leaf.key()).collect::<Vec<_>>();
    let mut nodes = pairing.sequences.nodes;
    complete_alignment_for_leaves(
      &nwk_parsed.graph,
      &mut nodes,
      4,
      &alphabet,
      MissingLeafLimit { ignored: false },
      "alignment 'aln.fasta'",
      &NoopProgress,
    )
    .unwrap();
    let leaf_seqs = leaves.iter().map(|key| nodes[key].seq.clone()).collect::<Vec<_>>();

    assert_eq!(
      (
        vec![
          Some(Seq::try_from_str("ACGT").unwrap()),
          Some(Seq::try_from_str("ACGT").unwrap()),
          Some(Seq::try_from_str("GGGG").unwrap()),
        ],
        vec![o!("A")],
      ),
      (leaf_seqs, pairing.duplicate_names)
    );
  }

  #[test]
  fn test_sanitize_to_alphabet_folds_stop_into_unknown_for_no_stop_alphabet() {
    let aa = Alphabet::new(AlphabetName::Aa).unwrap();
    let aa_no_stop = Alphabet::new(AlphabetName::AaNoStop).unwrap();
    let seq = Seq::try_from_str("MC*X-").unwrap();

    let (kept, changed_aa) = sanitize_to_alphabet(&seq, &aa);
    assert_eq!(0, changed_aa);
    assert_eq!(seq, kept);

    let (folded, changed) = sanitize_to_alphabet(&seq, &aa_no_stop);
    assert_eq!(1, changed);
    assert_eq!(Seq::try_from_str("MCXX-").unwrap(), folded);
  }

  mod helpers {
    use super::*;

    pub(super) fn two_leaf_tree() -> (Graph, NodeNames) {
      let nwk_parsed = nwk_read(b"(A:0.1,B:0.1)root;".as_slice()).unwrap();
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      (graph, names)
    }

    pub(super) fn three_leaf_tree() -> (Graph, NodeNames) {
      let nwk_parsed = nwk_read(b"(A:0.1,B:0.1,C:0.1)root;".as_slice()).unwrap();
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      (graph, names)
    }

    pub(super) fn four_leaf_tree() -> (Graph, NodeNames) {
      let nwk_parsed = nwk_read(b"((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1)root;".as_slice()).unwrap();
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      (graph, names)
    }

    pub(super) fn records(entries: &[(&str, &str)]) -> Vec<AlignmentRecord> {
      entries
        .iter()
        .map(|(name, seq)| AlignmentRecord {
          name: (*name).to_owned(),
          seq: Seq::try_from_str(seq).unwrap(),
        })
        .collect()
    }

    pub(super) fn complete(
      graph: &Graph,
      names: &NodeNames,
      records: Vec<AlignmentRecord>,
      alphabet: &Alphabet,
      limit: MissingLeafLimit,
    ) -> Result<BTreeMap<String, Seq>, Report> {
      let mut sequences = pair_leaf_sequences(graph, names, records).sequences;
      let length = sequences.common_length()?;
      complete_alignment_for_leaves(
        graph,
        &mut sequences.nodes,
        length,
        alphabet,
        limit,
        "alignment 'aln.fasta'",
        &NoopProgress,
      )?;
      Ok(
        sequences
          .nodes
          .into_values()
          .filter_map(|node| Some((node.name, node.seq?)))
          .collect(),
      )
    }
  }
}
