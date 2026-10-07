#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{F81Params, JC69Params, f81, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts};
  use crate::optimize::params::ExistingBranchLengths;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::{DenseReconstruction, MarginalReconstruction};
  use crate::test_utils::leaf_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use ndarray::array;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::{FastaRecord, fasta_read};
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn biased_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    fasta_read(
      indoc! {r#"
        >A
        AAAAAAAAAAAAAAAA
        >B
        CCCCCCCCCCCCCCCC
        >C
        AAAAAAAAAAAAAAAA
        >D
        CCCCCCCCCCCCCCCC
      "#}
      .as_bytes(),
      &alphabet,
    )
  }

  fn setup_dense_jc69(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[AlignmentRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<MarginalReconstruction, Report> {
    let alphabet = Alphabet::default();
    let partition = PartitionMarginalDense::new(alphabet, graph, &leaf_seq_inputs(graph, names, aln.to_vec()))?;
    let reconstruction =
      MarginalReconstruction::Dense(DenseReconstruction::seeded(partition, jc69(JC69Params::default())?));
    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;
    Ok(reconstruction)
  }

  fn get_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .map(|edge| branch_lengths[&edge.key()].unwrap_or(0.0))
      .collect()
  }

  #[test]
  fn test_stale_jc69_messages_bias_initial_guess() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = biased_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let f81_gtr = f81(F81Params {
      pi: Some(array![0.7, 0.1, 0.1, 0.1]),
      ..Default::default()
    })?;

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_stale_names = nwk_parsed.names();
    let graph_stale = nwk_parsed.graph;
    let mut branch_lengths_stale = nwk_parsed.branch_lengths;
    let reconstruction_stale = setup_dense_jc69(&graph_stale, &graph_stale_names, &aln, &branch_lengths_stale)?;
    let (mut reconstruction_stale, _) =
      reconstruction_stale.marginal_update(&graph_stale, &branch_lengths_or_zero(&branch_lengths_stale))?;
    *reconstruction_stale.gtr_mut() = f81_gtr.clone();
    {
      let total_length = reconstruction_stale.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph_stale, &reconstruction_stale);
      let sub_counts = gather_edge_sub_counts(&graph_stale, &reconstruction_stale)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_stale, &reconstruction_stale)?;
      initial_guess_mixed(
        &graph_stale,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths_stale,
      )?;
    }
    let bl_stale = get_branch_lengths(&graph_stale, &branch_lengths_stale);

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph_fresh_names = nwk_parsed.names();
    let graph_fresh = nwk_parsed.graph;
    let mut branch_lengths_fresh = nwk_parsed.branch_lengths;
    let reconstruction_fresh = setup_dense_jc69(&graph_fresh, &graph_fresh_names, &aln, &branch_lengths_fresh)?;
    let (mut reconstruction_fresh, _) =
      reconstruction_fresh.marginal_update(&graph_fresh, &branch_lengths_or_zero(&branch_lengths_fresh))?;
    *reconstruction_fresh.gtr_mut() = f81_gtr;
    let (reconstruction_fresh, _) =
      reconstruction_fresh.marginal_update(&graph_fresh, &branch_lengths_or_zero(&branch_lengths_fresh))?;
    {
      let total_length = reconstruction_fresh.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph_fresh, &reconstruction_fresh);
      let sub_counts = gather_edge_sub_counts(&graph_fresh, &reconstruction_fresh)?;
      let effective_lengths = gather_edge_effective_lengths(&graph_fresh, &reconstruction_fresh)?;
      initial_guess_mixed(
        &graph_fresh,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths_fresh,
      )?;
    }
    let bl_fresh = get_branch_lengths(&graph_fresh, &branch_lengths_fresh);

    assert_ne!(
      bl_stale, bl_fresh,
      "Stale JC69 messages should produce different initial branch lengths than fresh F81 messages"
    );

    Ok(())
  }

  #[test]
  fn test_initial_guess_idempotent_after_gtr_update() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = biased_alignment()?.into_iter().map(AlignmentRecord::from).collect();
    let f81_gtr = f81(F81Params {
      pi: Some(array![0.7, 0.1, 0.1, 0.1]),
      ..Default::default()
    })?;

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_dense_jc69(&graph, &names, &aln, &branch_lengths)?;
    let (mut reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    *reconstruction.gtr_mut() = f81_gtr.clone();
    let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    {
      let total_length = reconstruction.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
      let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
      initial_guess_mixed(
        &graph,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths,
      )?;
    }
    let bl_first = get_branch_lengths(&graph, &branch_lengths);

    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let graph2_names = nwk_parsed.names();
    let graph2 = nwk_parsed.graph;
    let mut branch_lengths2 = nwk_parsed.branch_lengths;
    let reconstruction2 = setup_dense_jc69(&graph2, &graph2_names, &aln, &branch_lengths2)?;
    let (mut reconstruction2, _) =
      reconstruction2.marginal_update(&graph2, &branch_lengths_or_zero(&branch_lengths2))?;
    *reconstruction2.gtr_mut() = f81_gtr;
    let (reconstruction2, _) = reconstruction2.marginal_update(&graph2, &branch_lengths_or_zero(&branch_lengths2))?;
    {
      let total_length = reconstruction2.sequence_length();
      let indel_counts = gather_edge_indel_counts(&graph2, &reconstruction2);
      let sub_counts = gather_edge_sub_counts(&graph2, &reconstruction2)?;
      let effective_lengths = gather_edge_effective_lengths(&graph2, &reconstruction2)?;
      initial_guess_mixed(
        &graph2,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        ExistingBranchLengths::Overwrite,
        false,
        &mut branch_lengths2,
      )?;
    }
    let bl_second = get_branch_lengths(&graph2, &branch_lengths2);

    assert_eq!(
      bl_first, bl_second,
      "Same initialization sequence should produce identical branch lengths"
    );

    Ok(())
  }
}
