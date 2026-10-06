#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::GtrModelName;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::create::Representation;
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::partition::marginal::reconstruction::SparseReconstruction;
  use crate::partition::marginal::sample::SampleMode;
  use crate::seq::mutation::{MutationTrack, stream_sequence_mutations};
  use crate::test_utils::RecordingSeqSink;
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_utils::sync::random::get_random_number_generator;

  #[test]
  fn test_partition_timetree_edge_sub_count_unavailable_before_inference() -> Result<(), Report> {
    let (graph, seeded) = helpers::seeded_sparse()?;
    let partition = MarginalReconstruction::Sparse(seeded);
    for edge in graph.get_edges() {
      assert_eq!(None, partition.edge_sub_count(&graph, edge.key())?);
    }
    Ok(())
  }

  #[test]
  fn test_partition_timetree_edge_sub_count_matches_edge_subs_after_inference() -> Result<(), Report> {
    let (graph, seeded) = helpers::seeded_sparse()?;
    let (updated, _) = seeded.marginal_update(&graph, &branch_lengths_or_zero(&helpers::branch_lengths()?))?;
    let partition = MarginalReconstruction::Sparse(updated);
    for edge in graph.get_edges() {
      let expected = partition.edge_subs(&graph, edge.key())?.len();
      assert_eq!(Some(expected), partition.edge_sub_count(&graph, edge.key())?);
    }
    Ok(())
  }

  #[rstest]
  #[case::dense(Representation::Dense)]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_sample_sequences_argmax_samples_no_node(#[case] representation: Representation) -> Result<(), Report> {
    let (graph, names, reconstruction) = helpers::updated_four_leaves(representation)?;
    let keys = helpers::keys(&graph, &names, &["root", "X", "Y"]);

    let actual = reconstruction.sample_sequences(&graph, SampleMode::Argmax, &mut get_random_number_generator(0))?;

    assert_eq!(BTreeMap::new(), actual);
    for key in keys {
      assert_eq!("ACGT", reconstruction.node_sequence(&graph, false, key)?.to_string());
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::internal_only(false, &[
    ("root", true), ("X", true), ("A", false), ("B", false), ("Y", true), ("C", false), ("D", false),
  ])]
  #[case::with_leaves(true, &[
    ("root", true), ("X", true), ("A", true), ("B", true), ("Y", true), ("C", true), ("D", true),
  ])]
  #[trace]
  fn test_stream_sequences_visits_all_nodes_in_preorder_and_flags_emitted_ones(
    #[values(Representation::Dense, Representation::Sparse)] representation: Representation,
    #[case] include_leaves: bool,
    #[case] expected: &[(&str, bool)],
  ) -> Result<(), Report> {
    let (graph, names, reconstruction) = helpers::updated_four_leaves(representation)?;
    let mut sink = RecordingSeqSink::default();

    stream_sequence_mutations(
      &graph,
      reconstruction.alphabet(),
      &MutationTrack::Nucleotide,
      include_leaves,
      true,
      |key| reconstruction.node_sequence(&graph, false, key),
      |edge_key| reconstruction.edge_indels(edge_key),
      Some(&mut sink),
    )?;

    let expected = expected
      .iter()
      .map(|&(name, emitted)| (name.to_owned(), emitted))
      .collect::<Vec<_>>();
    let actual = sink
      .items
      .into_iter()
      .map(|(key, emitted, _)| (names[&key].clone().unwrap(), emitted))
      .collect::<Vec<_>>();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dense_root( Representation::Dense,  SampleMode::Root, &["root"])]
  #[case::dense_all(  Representation::Dense,  SampleMode::All,  &["root", "X", "Y"])]
  #[case::sparse_root(Representation::Sparse, SampleMode::Root, &["root"])]
  #[case::sparse_all( Representation::Sparse, SampleMode::All,  &["root", "X", "Y"])]
  #[trace]
  fn test_reconstruct_sequences_sample_mode_selects_sampled_internal_nodes(
    #[case] representation: Representation,
    #[case] mode: SampleMode,
    #[case] expected_sampled: &[&str],
  ) -> Result<(), Report> {
    let (graph, names, reconstruction) = helpers::updated_four_leaves(representation)?;

    let actual = reconstruction.sample_sequences(&graph, mode, &mut get_random_number_generator(0))?;

    let mut expected_keys = helpers::keys(&graph, &names, expected_sampled);
    expected_keys.sort_unstable();
    assert_eq!(expected_keys, actual.keys().copied().collect::<Vec<_>>());
    for seq in actual.values() {
      assert_eq!(4, seq.len());
      assert!(
        seq.iter().all(|state| b"ACGT".contains(&u8::from(*state))),
        "a sampled state must lie in the support of the posterior profile, found {seq}"
      );
      if matches!(representation, Representation::Sparse) {
        assert!(seq.as_str().starts_with("ACG"), "sparse sampling keeps invariant columns, found {seq}");
      }
    }
    Ok(())
  }

  mod helpers {
    use super::*;
    use crate::partition::create::build_marginal_partition;
    use crate::progress::NoopProgress;
    use crate::test_utils::find_node_key_by_name;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::node::GraphNodeKey;

    const TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    pub(super) fn seeded_sparse() -> Result<(Graph, SparseReconstruction), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = fasta_read(
        indoc! {br#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#}
        .as_slice(),
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let nwk_parsed = nwk_read(TREE.as_bytes())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let fitch = create_fitch_partition(&graph, 0, alphabet, leaf_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let gtr = jc69(JC69Params::default())?;
      Ok((graph, SparseReconstruction::seeded(partition, gtr, node_states)))
    }

    pub(super) fn branch_lengths() -> Result<BTreeMap<GraphEdgeKey, Option<f64>>, Report> {
      Ok(nwk_read(TREE.as_bytes())?.branch_lengths)
    }

    pub(super) fn updated_four_leaves(
      representation: Representation,
    ) -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>, MarginalReconstruction), Report> {
      let nwk_parsed = nwk_read(b"((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)Y:0.1)root;".as_slice())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let aln: Vec<AlignmentRecord> = fasta_read(
        indoc! {br#"
        >A
        ACGT
        >B
        ACGT
        >C
        ACGT
        >D
        ACGA
      "#}
        .as_slice(),
        &Alphabet::new(AlphabetName::Nuc)?,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let reconstruction = build_marginal_partition(
        representation,
        GtrModelName::JC69,
        &graph,
        0,
        Alphabet::new(AlphabetName::Nuc)?,
        leaf_seq_inputs(&graph, &names, aln),
        &branch_lengths,
        &NoopProgress,
      )?;
      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths)?;
      Ok((graph, names, reconstruction))
    }

    pub(super) fn keys(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      node_names: &[&str],
    ) -> Vec<GraphNodeKey> {
      node_names
        .iter()
        .map(|name| find_node_key_by_name(graph, names, name).expect("node must exist"))
        .collect()
    }
  }
}
