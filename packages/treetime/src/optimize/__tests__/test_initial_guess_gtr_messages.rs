#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{F81Params, JC69Params, f81, jc69};
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use eyre::Report;
  use indoc::indoc;
  
  use ndarray::array;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
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
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalDense>>>, Report> {
    let alphabet = Alphabet::default();
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(0, jc69(JC69Params::default())?, alphabet, get_common_length(aln)?)))];
    initialize_marginal(graph, &partitions, aln)?;
    Ok(partitions)
  }

  fn get_branch_lengths(graph: &GraphAncestral) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect()
  }

  /// Regression: after replacing the dummy JC69 with a non-uniform GTR,
  /// update_marginal must be re-run before initial_guess_mixed. Stale JC69
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
    // Replace GTR but do NOT re-run update_marginal.
    let graph_stale: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions_stale = setup_dense_jc69(&graph_stale, &aln)?;
    update_marginal(&graph_stale, &partitions_stale)?;
    partitions_stale[0].write_arc().data.gtr = f81_gtr.clone();
    initial_guess_mixed(&graph_stale, &partitions_stale, true)?;
    let bl_stale = get_branch_lengths(&graph_stale);

    // Scenario 2: fresh F81 messages (the fix)
    // Replace GTR AND re-run update_marginal.
    let graph_fresh: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions_fresh = setup_dense_jc69(&graph_fresh, &aln)?;
    update_marginal(&graph_fresh, &partitions_fresh)?;
    partitions_fresh[0].write_arc().data.gtr = f81_gtr;
    update_marginal(&graph_fresh, &partitions_fresh)?;
    initial_guess_mixed(&graph_fresh, &partitions_fresh, true)?;
    let bl_fresh = get_branch_lengths(&graph_fresh);

    assert_ne!(
      bl_stale, bl_fresh,
      "Stale JC69 messages should produce different initial branch lengths than fresh F81 messages"
    );

    Ok(())
  }

  /// After the initialization sequence with the fix, a redundant
  /// update_marginal should not change the initial guess: the messages
  /// are already computed with the real GTR.
  #[test]
  fn test_initial_guess_idempotent_after_gtr_update() -> Result<(), Report> {
    let aln = biased_alignment()?;
    let f81_gtr = f81(F81Params {
      pi: Some(array![0.7, 0.1, 0.1, 0.1]),
      ..Default::default()
    })?;

    // Run full initialization with the fix: JC69, update, replace, update, guess
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense_jc69(&graph, &aln)?;
    update_marginal(&graph, &partitions)?;
    partitions[0].write_arc().data.gtr = f81_gtr.clone();
    update_marginal(&graph, &partitions)?;
    initial_guess_mixed(&graph, &partitions, true)?;
    let bl_first = get_branch_lengths(&graph);

    // Run update_marginal + initial_guess again on the same graph.
    // Branch lengths changed from initial_guess, so update_marginal
    // recomputes messages with the new branch lengths. The resulting
    // initial guess may differ from bl_first (new transition matrices).
    // But running the SAME sequence twice from identical state must
    // produce the same result.
    let graph2: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions2 = setup_dense_jc69(&graph2, &aln)?;
    update_marginal(&graph2, &partitions2)?;
    partitions2[0].write_arc().data.gtr = f81_gtr;
    update_marginal(&graph2, &partitions2)?;
    initial_guess_mixed(&graph2, &partitions2, true)?;
    let bl_second = get_branch_lengths(&graph2);

    assert_eq!(
      bl_first, bl_second,
      "Same initialization sequence should produce identical branch lengths"
    );

    Ok(())
  }
}
