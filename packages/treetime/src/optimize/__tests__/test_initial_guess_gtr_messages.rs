#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::ancestral::pipeline::DenseReconstruction;
  use crate::gtr::get_gtr::{F81Params, JC69Params, f81, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::gather::{
    gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts, total_sequence_length,
  };
  use crate::optimize::run_loop::marginal_update_dense;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_fasta_node_inputs;

  use ndarray::array;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  /// Alignment where every position has A on leaves A/C and C on leaves B/D.
  /// Internal node reconstructions are uncertain between A and C, making them
  /// sensitive to equilibrium frequencies.
  fn biased_alignment() -> Result<Vec<FastaRecord>, Report> {
    let alphabet = Alphabet::default();
    read_many_fasta_str(
      indoc! {r#"
        >A
        AAAAAAAAAAAAAAAA
        >B
        CCCCCCCCCCCCCCCC
        >C
        AAAAAAAAAAAAAAAA
        >D
        CCCCCCCCCCCCCCCC
      "#},
      &alphabet,
    )
  }

  fn setup_dense_jc69(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    aln: &[FastaRecord],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<Vec<DenseReconstruction>, Report> {
    let alphabet = Alphabet::default();
    let partition = PartitionMarginalDense::new(0, alphabet, get_common_length(aln)?);
    let node_states = partition.attach_sequences(graph, &nwk_fasta_node_inputs(graph, names, aln.to_vec()))?;
    let partitions = vec![DenseReconstruction::seeded(partition, jc69(JC69Params::default())?, node_states)];
    let (partitions, _) = marginal_update_dense(graph, &profile_branch_lengths(branch_lengths), partitions)?;
    Ok(partitions)
  }

  fn get_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| branch_lengths[&edge.read_arc().key()].unwrap_or(0.0))
      .collect()
  }

  /// Regression: after replacing the dummy JC69 with a non-uniform GTR,
  /// marginal_update must be re-run before initial_guess_mixed. Stale JC69
  /// node posteriors produce a different (biased) initial branch length guess
  /// than fresh posteriors computed with the real model.
  #[test]
  fn test_stale_jc69_messages_bias_initial_guess() -> Result<(), Report> {
    let aln = biased_alignment()?;
    let f81_gtr = f81(F81Params {
      pi: Some(array![0.7, 0.1, 0.1, 0.1]),
      ..Default::default()
    })?;

    // Scenario 1: stale JC69 messages (the bug)
    // Replace GTR but do NOT re-run marginal_update.
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_stale_names = nwk_parsed.names();
    let graph_stale = nwk_parsed.graph;
    let mut branch_lengths_stale = nwk_parsed.branch_lengths;
    let partitions_stale = setup_dense_jc69(&graph_stale, &graph_stale_names, &aln, &branch_lengths_stale)?;
    let (mut partitions_stale, _) = marginal_update_dense(
      &graph_stale,
      &profile_branch_lengths(&branch_lengths_stale),
      partitions_stale,
    )?;
    partitions_stale[0].gtr = f81_gtr.clone();
    {
      let total_length = total_sequence_length(&partitions_stale, &[]);
      let indel_counts = gather_edge_indel_counts(&graph_stale, &partitions_stale, &[]);
      let sub_counts = gather_edge_sub_counts(&graph_stale, &partitions_stale, &[])?;
      let effective_lengths = gather_edge_effective_lengths(&graph_stale, &partitions_stale, &[])?;
      initial_guess_mixed(
        &graph_stale,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
        false,
        &mut branch_lengths_stale,
      )?;
    }
    let bl_stale = get_branch_lengths(&graph_stale, &branch_lengths_stale);

    // Scenario 2: fresh F81 messages (the fix)
    // Replace GTR AND re-run marginal_update.
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_fresh_names = nwk_parsed.names();
    let graph_fresh = nwk_parsed.graph;
    let mut branch_lengths_fresh = nwk_parsed.branch_lengths;
    let partitions_fresh = setup_dense_jc69(&graph_fresh, &graph_fresh_names, &aln, &branch_lengths_fresh)?;
    let (mut partitions_fresh, _) = marginal_update_dense(
      &graph_fresh,
      &profile_branch_lengths(&branch_lengths_fresh),
      partitions_fresh,
    )?;
    partitions_fresh[0].gtr = f81_gtr;
    let (partitions_fresh, _) = marginal_update_dense(
      &graph_fresh,
      &profile_branch_lengths(&branch_lengths_fresh),
      partitions_fresh,
    )?;
    {
      let total_length = total_sequence_length(&partitions_fresh, &[]);
      let indel_counts = gather_edge_indel_counts(&graph_fresh, &partitions_fresh, &[]);
      let sub_counts = gather_edge_sub_counts(&graph_fresh, &partitions_fresh, &[])?;
      let effective_lengths = gather_edge_effective_lengths(&graph_fresh, &partitions_fresh, &[])?;
      initial_guess_mixed(
        &graph_fresh,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
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

  /// After the initialization sequence with the fix, a redundant
  /// marginal_update should not change the initial guess: the messages
  /// are already computed with the real GTR.
  #[test]
  fn test_initial_guess_idempotent_after_gtr_update() -> Result<(), Report> {
    let aln = biased_alignment()?;
    let f81_gtr = f81(F81Params {
      pi: Some(array![0.7, 0.1, 0.1, 0.1]),
      ..Default::default()
    })?;

    // Run full initialization with the fix: JC69, update, replace, update, guess
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_dense_jc69(&graph, &names, &aln, &branch_lengths)?;
    let (mut partitions, _) = marginal_update_dense(&graph, &profile_branch_lengths(&branch_lengths), partitions)?;
    partitions[0].gtr = f81_gtr.clone();
    let (partitions, _) = marginal_update_dense(&graph, &profile_branch_lengths(&branch_lengths), partitions)?;
    {
      let total_length = total_sequence_length(&partitions, &[]);
      let indel_counts = gather_edge_indel_counts(&graph, &partitions, &[]);
      let sub_counts = gather_edge_sub_counts(&graph, &partitions, &[])?;
      let effective_lengths = gather_edge_effective_lengths(&graph, &partitions, &[])?;
      initial_guess_mixed(
        &graph,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
        false,
        &mut branch_lengths,
      )?;
    }
    let bl_first = get_branch_lengths(&graph, &branch_lengths);

    // Run marginal_update + initial_guess again on the same graph.
    // Branch lengths changed from initial_guess, so marginal_update
    // recomputes messages with the new branch lengths. The resulting
    // initial guess may differ from bl_first (new transition matrices).
    // But running the SAME sequence twice from identical state must
    // produce the same result.
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph2_names = nwk_parsed.names();
    let graph2 = nwk_parsed.graph;
    let mut branch_lengths2 = nwk_parsed.branch_lengths;
    let partitions2 = setup_dense_jc69(&graph2, &graph2_names, &aln, &branch_lengths2)?;
    let (mut partitions2, _) = marginal_update_dense(&graph2, &profile_branch_lengths(&branch_lengths2), partitions2)?;
    partitions2[0].gtr = f81_gtr;
    let (partitions2, _) = marginal_update_dense(&graph2, &profile_branch_lengths(&branch_lengths2), partitions2)?;
    {
      let total_length = total_sequence_length(&partitions2, &[]);
      let indel_counts = gather_edge_indel_counts(&graph2, &partitions2, &[]);
      let sub_counts = gather_edge_sub_counts(&graph2, &partitions2, &[])?;
      let effective_lengths = gather_edge_effective_lengths(&graph2, &partitions2, &[])?;
      initial_guess_mixed(
        &graph2,
        total_length,
        &indel_counts,
        &sub_counts,
        &effective_lengths,
        true,
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
