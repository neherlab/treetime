#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::optimize_unified::{initial_guess_mixed, run_optimize_mixed};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use std::sync::LazyLock;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

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

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    let max_iter = 50;

    let mut lh_history = Vec::with_capacity(max_iter);
    lh_history.push(initial_lh);

    for _ in 0..max_iter {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      let lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
      lh_history.push(lh);
    }

    let final_lh = lh_history[lh_history.len() - 1];

    // Verify final likelihood is in expected range for this tree/alignment
    // JC69 on 16-site alignment with 4 leaves
    assert!(
      final_lh > -100.0 && final_lh < -30.0,
      "Final log-lh {final_lh} outside expected range [-100, -30]"
    );

    // Check convergence: variance over last 5 iterations should be small
    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean = last_5.iter().sum::<f64>() / 5.0;
    let variance = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 5.0;
    assert!(
      variance < 1.0,
      "Optimization should stabilize: variance of last 5 iterations = {variance}"
    );

    Ok(())
  }

  #[test]
  fn test_optimization_improves_or_maintains_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    // Initial log-lh should be negative (log of probability < 1)
    assert!(initial_lh < 0.0, "Initial log-LH should be negative: {initial_lh}");

    // Run several optimization steps
    for _ in 0..10 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final likelihood should be at least as good as initial (within small tolerance)
    // Note: individual iterations may fluctuate, but overall trend should improve or maintain
    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization significantly worsened likelihood: {initial_lh} -> {final_lh}"
    );

    // Final log-lh should be in reasonable range for JC69 on this alignment
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH unreasonably low: {final_lh}");

    Ok(())
  }

  #[test]
  fn test_optimization_produces_valid_branch_lengths() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Collect initial branch lengths
    let initial_total: f64 = graph
      .get_edges()
      .iter()
      .filter_map(|e| e.read_arc().payload().read_arc().branch_length())
      .sum();

    // Run several optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    // Verify all branch lengths are in valid range
    let mut final_total = 0.0;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let branch_length = edge.payload().read_arc().branch_length();
      if let Some(bl) = branch_length {
        // Branch lengths should be non-negative and reasonable (< 10 subs/site)
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
        final_total += bl;
      }
    }

    // Total tree length should be reasonable (not zero, not huge)
    assert!(final_total > 0.0, "Total tree length should be positive");
    assert!(
      final_total < 50.0,
      "Total tree length unreasonably large: {final_total}"
    );

    // Tree length should be in same order of magnitude as initial
    // (optimization shouldn't drastically change overall scale)
    assert!(
      final_total > initial_total * 0.1 && final_total < initial_total * 10.0,
      "Tree length changed too drastically: {initial_total} -> {final_total}"
    );

    Ok(())
  }

  // ==========================================================================
  // Idempotence tests
  // ==========================================================================

  #[test]
  fn test_optimization_converges_with_valid_branch_lengths() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    let mut lh_history = Vec::with_capacity(20);

    // Run optimization iterations
    for i in 0..20 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      let lh = update_marginal(&graph, &dense_partitions)? + update_marginal(&graph, &sparse_partitions)?;

      lh_history.push(lh);

      // After each iteration, all branch lengths should be non-negative and bounded
      for edge in graph.get_edges() {
        let branch_length = edge.read_arc().payload().read_arc().branch_length();
        if let Some(bl) = branch_length {
          assert!(bl >= 0.0, "Branch length should be non-negative at iter {i}: {bl}");
          assert!(bl < 10.0, "Branch length too large at iter {i}: {bl}");
        }
      }
    }

    // Check convergence: variance over last 5 iterations should be small
    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean = last_5.iter().sum::<f64>() / 5.0;
    let variance = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 5.0;
    assert!(
      variance < 1.0,
      "Optimization should stabilize: variance of last 5 iterations = {variance}"
    );

    // Final likelihood should be in expected range
    let final_lh = lh_history[19];
    assert!(final_lh < 0.0, "Final log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Final log-LH unreasonably low: {final_lh}");

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
    // Tree with zero branch length on edge to A
    let tree_newick = "((A:0.0,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Run multiple optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final log-lh should be in reasonable range
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH should be reasonable: {final_lh}");

    // Branch lengths should be non-negative and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }

  #[test]
  fn test_optimization_handles_very_short_branches() -> Result<(), Report> {
    // Tree with very short branch lengths (all 0.0001)
    let tree_newick = "((A:0.0001,B:0.0001)AB:0.0001,(C:0.0001,D:0.0001)CD:0.0001)root:0.0001;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Run optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final log-lh should be negative and reasonable
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH should be reasonable: {final_lh}");

    // Branch lengths should be non-negative and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 10.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }

  #[test]
  fn test_optimization_handles_long_branches() -> Result<(), Report> {
    // Tree with longer branch lengths (some > 1 sub/site)
    let tree_newick = "((A:1.0,B:2.0)AB:1.0,(C:2.0,D:1.2)CD:0.5)root:0.1;";
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(tree_newick)?;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &aln)?;

    // Run optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &dense_partitions, &sparse_partitions)?;
      update_marginal(&graph, &dense_partitions)?;
      update_marginal(&graph, &sparse_partitions)?;
    }

    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Final log-lh should be negative and reasonable
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -200.0, "Log-LH should be reasonable: {final_lh}");

    // Branch lengths should be non-negative and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl >= 0.0, "Branch length should be non-negative: {bl}");
        assert!(bl < 20.0, "Branch length unreasonably large: {bl}");
      }
    }

    Ok(())
  }
}
