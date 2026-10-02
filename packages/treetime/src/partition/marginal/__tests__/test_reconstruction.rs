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
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::mutation::{MutationTrack, stream_sequence_mutations};
  use crate::test_utils::RecordingSeqSink;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rand::SeedableRng;
  use rand::rngs::StdRng;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

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

    let actual = reconstruction.sample_sequences(&graph, SampleMode::Argmax, &mut StdRng::seed_from_u64(0))?;

    assert_eq!(BTreeMap::new(), actual);
    for key in keys {
      assert_eq!("ACGT", reconstruction.node_sequence(&graph, false, key)?.to_string());
    }
    Ok(())
  }

  #[rstest]
  #[case::dense(Representation::Dense)]
  #[case::sparse(Representation::Sparse)]
  #[trace]
  fn test_stream_sequences_visits_all_nodes_in_preorder_and_flags_emitted_ones(
    #[case] representation: Representation,
    #[values(false, true)] include_leaves: bool,
  ) -> Result<(), Report> {
    let (graph, names, reconstruction) = helpers::updated_four_leaves(representation)?;
    let mut sink = RecordingSeqSink::default();

    stream_sequence_mutations(
      &graph,
      reconstruction.alphabet(),
      &MutationTrack::Nucleotide,
      include_leaves,
      |key| reconstruction.node_sequence(&graph, false, key),
      |edge_key| reconstruction.edge_indels(edge_key),
      Some(&mut sink),
    )?;

    let expected = ["root", "X", "A", "B", "Y", "C", "D"]
      .into_iter()
      .map(|name| {
        (
          name.to_owned(),
          include_leaves || !matches!(name, "A" | "B" | "C" | "D"),
        )
      })
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

    let actual = reconstruction.sample_sequences(&graph, mode, &mut StdRng::seed_from_u64(0))?;

    let mut expected_keys = helpers::keys(&graph, &names, expected_sampled);
    expected_keys.sort_unstable();
    assert_eq!(expected_keys, actual.keys().copied().collect::<Vec<_>>());
    for seq in actual.values() {
      assert_eq!(4, seq.len());
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
      let aln: Vec<AlignmentRecord> = read_many_fasta_str(
        indoc! {r#"
        >A
        ACATCGCCNNA--GAC
        >B
        GCATCCCTGTA-NG--
        >C
        CCGGCGATGTRTTG--
        >D
        TCGGCCGTGTRTTG--
      "#},
        &alphabet,
      )?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();
      let nwk_parsed = nwk_read_str(TREE)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let fitch = create_fitch_partition(&graph, 0, alphabet, &node_seq_inputs(&graph, &names, aln))?;
      let (partition, node_states) = fitch.into_marginal_sparse(&graph)?;
      let gtr = jc69(JC69Params::default())?;
      Ok((graph, SparseReconstruction::seeded(partition, gtr, node_states)))
    }

    pub(super) fn branch_lengths() -> Result<BTreeMap<GraphEdgeKey, Option<f64>>, Report> {
      Ok(nwk_read_str(TREE)?.branch_lengths)
    }

    pub(super) fn updated_four_leaves(
      representation: Representation,
    ) -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>, MarginalReconstruction), Report> {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)X:0.1,(C:0.1,D:0.1)Y:0.1)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let aln: Vec<AlignmentRecord> = read_many_fasta_str(
        indoc! {r#"
        >A
        ACGT
        >B
        ACGT
        >C
        ACGT
        >D
        ACGA
      "#},
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
        &node_seq_inputs(&graph, &names, aln),
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
