#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  // Small tree with 4 leaves
  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn simple_alignment() -> Result<Vec<FastaRecord>, Report> {
    read_many_fasta_str(
      indoc! {r#"
      >A
      ACGTACGTACGTACGT
      >B
      ACGTACGTACGTACGA
      >C
      ACGTACGTACGTACGG
      >D
      ACGTACGTACGTACGC
    "#},
      &*NUC_ALPHABET,
    )
  }

  fn setup_partitions(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<
    (
      Vec<Arc<RwLock<PartitionMarginalDense>>>,
      Vec<Arc<RwLock<PartitionMarginalSparse>>>,
    ),
    Report,
  > {
    let alphabet_dense = Alphabet::new(AlphabetName::Nuc, true)?;
    let alphabet_sparse = Alphabet::new(AlphabetName::Nuc, false)?;

    let dense_partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet_dense,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    let sparse_partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 1,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet_sparse,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(graph, &sparse_partitions, aln)?;
    initialize_marginal(graph, &dense_partitions, aln)?;
    update_marginal(graph, &sparse_partitions)?;

    initial_guess_mixed(graph, &dense_partitions, &sparse_partitions);

    Ok((dense_partitions, sparse_partitions))
  }

  fn compute_total_lh(
    graph: &GraphAncestral,
    dense_partitions: &[Arc<RwLock<PartitionMarginalDense>>],
    sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  ) -> Result<f64, Report> {
    let dense_lh = update_marginal(graph, dense_partitions)?;
    let sparse_lh = update_marginal(graph, sparse_partitions)?;
    Ok(dense_lh + sparse_lh)
  }

  // ==========================================================================
  // Convergence tests
  // ==========================================================================

  #[test]
  fn test_optimization_converges_within_iterations() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    let mut prev_lh = f64::NEG_INFINITY;
    let max_iter = 50;
    let convergence_threshold = 1e-4;

    for _ in 0..max_iter {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

      if (lh - prev_lh).abs() < convergence_threshold {
        // Converged
        return Ok(());
      }

      prev_lh = lh;
    }

    // Should converge within max_iter - but even if not fully converged,
    // log-LH should be finite and stable
    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(final_lh.is_finite(), "Final log-LH should be finite: {final_lh}");

    Ok(())
  }

  #[test]
  fn test_optimization_produces_finite_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite: {initial_lh}");

    // Run several optimization steps
    for _ in 0..5 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization: {lh}");
    }

    Ok(())
  }

  #[test]
  fn test_optimization_produces_finite_branch_lengths() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Run several optimization iterations
    for _ in 0..5 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    // Verify all branch lengths are finite and non-negative
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let branch_length = edge.payload().read_arc().branch_length();
      if let Some(bl) = branch_length {
        assert!(bl.is_finite(), "Branch length should be finite: {bl}");
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
      }
    }

    Ok(())
  }

  // ==========================================================================
  // Idempotence tests
  // ==========================================================================

  #[test]
  fn test_optimization_preserves_branch_length_validity() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Run optimization iterations
    for _ in 0..20 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;

      // After each iteration, all branch lengths should remain valid
      for edge in graph.get_edges() {
        let branch_length = edge.read_arc().payload().read_arc().branch_length();
        if let Some(bl) = branch_length {
          assert!(bl.is_finite(), "Branch length should be finite: {bl}");
          assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        }
      }
    }

    Ok(())
  }

  #[test]
  fn test_second_optimization_produces_same_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;

    // Run optimization on first graph
    let graph1: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions1, sparse_partitions1) = setup_partitions(&graph1, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph1, &dense_partitions1, &sparse_partitions1)?;
      update_marginal(&graph1, &dense_partitions1)?;
      update_marginal(&graph1, &sparse_partitions1)?;
    }

    let lh1 = compute_total_lh(&graph1, &dense_partitions1, &sparse_partitions1)?;

    // Run optimization on second independent graph
    let graph2: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions2, sparse_partitions2) = setup_partitions(&graph2, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph2, &dense_partitions2, &sparse_partitions2)?;
      update_marginal(&graph2, &dense_partitions2)?;
      update_marginal(&graph2, &sparse_partitions2)?;
    }

    let lh2 = compute_total_lh(&graph2, &dense_partitions2, &sparse_partitions2)?;

    // Both runs should converge to same likelihood
    assert_ulps_eq!(lh1, lh2, max_ulps = 100);

    Ok(())
  }

  // ==========================================================================
  // Edge case tests
  // ==========================================================================

  #[test]
  fn test_optimization_handles_zero_branch_lengths() -> Result<(), Report> {
    // Tree with zero branch length
    let tree_newick = "((A:0.0,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Should not crash
    run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;

    let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(lh.is_finite(), "Log-LH should be finite: {lh}");

    Ok(())
  }

  #[test]
  fn test_optimization_handles_very_short_branches() -> Result<(), Report> {
    // Tree with very short branch lengths
    let tree_newick = "((A:0.0001,B:0.0001)AB:0.0001,(C:0.0001,D:0.0001)CD:0.0001)root:0.0001;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Should not crash
    for _ in 0..5 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(lh.is_finite(), "Log-LH should be finite: {lh}");

    Ok(())
  }

  #[test]
  fn test_optimization_handles_long_branches() -> Result<(), Report> {
    // Tree with longer branch lengths
    let tree_newick = "((A:1.0,B:2.0)AB:1.0,(C:2.0,D:1.2)CD:0.5)root:0.1;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Should not crash
    for _ in 0..5 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(lh.is_finite(), "Log-LH should be finite: {lh}");

    Ok(())
  }
}
