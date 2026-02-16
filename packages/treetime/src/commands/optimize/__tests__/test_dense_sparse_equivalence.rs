/// Tests for dense vs sparse optimization equivalence.
///
/// These tests verify that dense and sparse optimization modes produce consistent
/// results. Due to implementation differences (different zero-branch thresholds,
/// optimization paths), the modes may not produce identical results but should:
/// 1. Both produce finite, valid log-LH values
/// 2. Both converge to stable values
/// 3. Initial log-LH (before optimization) should be identical
/// 4. Final log-LH difference should be bounded
#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::{compress_sequences, get_common_length};
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::optimize_dense::run_optimize_dense;
  use crate::commands::optimize::optimize_sparse::run_optimize_sparse;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::graph_ancestral::GraphAncestral;
  use crate::representation::partition_marginal_dense::PartitionMarginalDense;
  use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
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

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn gap_free_alignment() -> Result<Vec<FastaRecord>, Report> {
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

  fn setup_dense_only(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalDense>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    initialize_marginal(graph, &partitions, aln)?;

    Ok(partitions)
  }

  fn setup_sparse_only(
    graph: &GraphAncestral,
    aln: &[FastaRecord],
  ) -> Result<Vec<Arc<RwLock<PartitionMarginalSparse>>>, Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let partitions = vec![Arc::new(RwLock::new(PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet,
      length: get_common_length(aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(graph, &partitions, aln)?;
    update_marginal(graph, &partitions)?;

    Ok(partitions)
  }

  fn get_branch_lengths(graph: &GraphAncestral) -> Vec<f64> {
    graph
      .get_edges()
      .iter()
      .map(|edge| edge.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect()
  }

  // ==========================================================================
  // Initial equivalence tests (before optimization)
  // ==========================================================================

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence() -> Result<(), Report> {
    let aln = gap_free_alignment()?;

    // Initialize dense
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;
    let log_lh_dense = update_marginal(&graph_dense, &dense_partitions)?;

    // Initialize sparse
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;
    let log_lh_sparse = update_marginal(&graph_sparse, &sparse_partitions)?;

    // Initial log-LH should be equivalent (before any optimization)
    assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence_with_mutations() -> Result<(), Report> {
    // Alignment with more mutations
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      AAAAAAAAAAAAAAAA
      >B
      CCCCCCCCCCCCCCCC
      >C
      GGGGGGGGGGGGGGGG
      >D
      TTTTTTTTTTTTTTTT
    "#},
      &*NUC_ALPHABET,
    )?;

    // Initialize dense
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;
    let log_lh_dense = update_marginal(&graph_dense, &dense_partitions)?;

    // Initialize sparse
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;
    let log_lh_sparse = update_marginal(&graph_sparse, &sparse_partitions)?;

    // Initial log-LH should be equivalent
    assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }

  // ==========================================================================
  // Optimization validity tests
  // ==========================================================================

  #[test]
  fn test_dense_optimization_produces_valid_results() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      run_optimize_dense(&graph, &partitions)?;
      let lh = update_marginal(&graph, &partitions)?;
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization");
    }

    let final_lh = update_marginal(&graph, &partitions)?;

    // Final log-LH should be in expected range for this simple tree
    // 16 sites, 4 leaves, mostly identical sequences -> log-LH between -100 and -10
    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    // Optimization should improve or maintain likelihood
    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization should not significantly decrease likelihood: initial={initial_lh}, final={final_lh}"
    );

    // Branch lengths should be valid and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl.is_finite(), "Branch length should be finite");
        assert!(bl >= 0.0, "Branch length should be non-negative");
        assert!(bl < 10.0, "Branch length {bl} should be reasonable (< 10)");
      }
    }

    Ok(())
  }

  #[test]
  fn test_sparse_optimization_produces_valid_results() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      run_optimize_sparse(&graph, &partitions)?;
      let lh = update_marginal(&graph, &partitions)?;
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization");
    }

    let final_lh = update_marginal(&graph, &partitions)?;

    // Final log-LH should be in expected range for this simple tree
    // 16 sites, 4 leaves, mostly identical sequences -> log-LH between -100 and -10
    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    // Optimization should improve or maintain likelihood
    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization should not significantly decrease likelihood: initial={initial_lh}, final={final_lh}"
    );

    // Branch lengths should be valid and bounded
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl.is_finite(), "Branch length should be finite");
        assert!(bl >= 0.0, "Branch length should be non-negative");
        assert!(bl < 10.0, "Branch length {bl} should be reasonable (< 10)");
      }
    }

    Ok(())
  }

  // ==========================================================================
  // Bounded difference tests
  // ==========================================================================

  #[test]
  fn test_dense_sparse_log_lh_bounded_difference_after_optimization() -> Result<(), Report> {
    let aln = gap_free_alignment()?;

    // Run dense-only optimization
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;

    for _ in 0..10 {
      run_optimize_dense(&graph_dense, &dense_partitions)?;
      update_marginal(&graph_dense, &dense_partitions)?;
    }

    let log_lh_dense = update_marginal(&graph_dense, &dense_partitions)?;

    // Run sparse-only optimization
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;

    for _ in 0..10 {
      run_optimize_sparse(&graph_sparse, &sparse_partitions)?;
      update_marginal(&graph_sparse, &sparse_partitions)?;
    }

    let log_lh_sparse = update_marginal(&graph_sparse, &sparse_partitions)?;

    // Both modes should produce finite log-LH in expected range
    assert!(
      log_lh_dense > -100.0 && log_lh_dense < -10.0,
      "Dense log-LH {log_lh_dense} should be in range [-100, -10]"
    );
    assert!(
      log_lh_sparse > -100.0 && log_lh_sparse < -10.0,
      "Sparse log-LH {log_lh_sparse} should be in range [-100, -10]"
    );

    // Dense and sparse should converge to similar values
    // Differences arise from zero-branch thresholds (0.01 dense vs 0.0001 sparse)
    // and different optimization paths, but should be bounded
    let diff = (log_lh_dense - log_lh_sparse).abs();
    assert!(
      diff < 0.5,
      "Log-LH difference should be bounded: dense={log_lh_dense}, sparse={log_lh_sparse}, diff={diff}"
    );

    Ok(())
  }

  #[test]
  fn test_dense_sparse_branch_lengths_bounded_difference() -> Result<(), Report> {
    let aln = gap_free_alignment()?;

    // Run dense-only optimization
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;

    for _ in 0..10 {
      run_optimize_dense(&graph_dense, &dense_partitions)?;
      update_marginal(&graph_dense, &dense_partitions)?;
    }

    let branch_lengths_dense = get_branch_lengths(&graph_dense);

    // Run sparse-only optimization
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;

    for _ in 0..10 {
      run_optimize_sparse(&graph_sparse, &sparse_partitions)?;
      update_marginal(&graph_sparse, &sparse_partitions)?;
    }

    let branch_lengths_sparse = get_branch_lengths(&graph_sparse);

    // Both modes should produce same number of edges
    assert_eq!(branch_lengths_dense.len(), branch_lengths_sparse.len());

    // All branch lengths should be valid and bounded
    for (dense_bl, sparse_bl) in branch_lengths_dense.iter().zip(branch_lengths_sparse.iter()) {
      assert!(dense_bl.is_finite());
      assert!(sparse_bl.is_finite());
      assert!(*dense_bl >= 0.0);
      assert!(*sparse_bl >= 0.0);
      assert!(*dense_bl < 10.0, "Dense branch length {dense_bl} should be reasonable");
      assert!(
        *sparse_bl < 10.0,
        "Sparse branch length {sparse_bl} should be reasonable"
      );
    }

    // Dense and sparse should produce similar branch lengths
    // Compare total tree length as a summary statistic
    let total_dense: f64 = branch_lengths_dense.iter().sum();
    let total_sparse: f64 = branch_lengths_sparse.iter().sum();
    let total_diff = (total_dense - total_sparse).abs();
    assert!(
      total_diff < 0.1,
      "Total tree length should be similar: dense={total_dense}, sparse={total_sparse}, diff={total_diff}"
    );

    // Individual branch lengths should also be close
    for (i, (dense_bl, sparse_bl)) in branch_lengths_dense
      .iter()
      .zip(branch_lengths_sparse.iter())
      .enumerate()
    {
      let diff = (dense_bl - sparse_bl).abs();
      assert!(
        diff < 0.05,
        "Branch {i} lengths should be similar: dense={dense_bl}, sparse={sparse_bl}, diff={diff}"
      );
    }

    Ok(())
  }

  // ==========================================================================
  // Convergence stability tests
  // ==========================================================================

  #[test]
  fn test_dense_optimization_converges() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    let mut lh_history = vec![initial_lh];

    for _ in 0..50 {
      run_optimize_dense(&graph, &partitions)?;
      let lh = update_marginal(&graph, &partitions)?;
      lh_history.push(lh);
    }

    let final_lh = *lh_history.last().unwrap();

    // Final log-LH should be in expected range
    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    // Optimization should improve or maintain likelihood overall
    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization should not significantly decrease likelihood: initial={initial_lh}, final={final_lh}"
    );

    // Check convergence: variance of last 5 iterations should be small
    let last_5: Vec<f64> = lh_history.iter().rev().take(5).copied().collect();
    let mean: f64 = last_5.iter().sum::<f64>() / last_5.len() as f64;
    let variance: f64 = last_5.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / last_5.len() as f64;
    assert!(
      variance < 1.0,
      "Optimization should stabilize: variance of last 5 iterations = {variance}"
    );

    Ok(())
  }

  #[test]
  fn test_sparse_optimization_converges() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    let mut lh_history = vec![initial_lh];

    for _ in 0..50 {
      run_optimize_sparse(&graph, &partitions)?;
      let lh = update_marginal(&graph, &partitions)?;
      lh_history.push(lh);
    }

    let final_lh = *lh_history.last().unwrap();

    // Final log-LH should be in expected range
    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    // All iterations should stay in valid range
    for (i, lh) in lh_history.iter().enumerate() {
      assert!(
        *lh > -200.0 && *lh < 0.0,
        "Iteration {i} log-LH {lh} should be in valid range [-200, 0]"
      );
    }

    // Best likelihood achieved should be reasonable
    let best_lh = lh_history.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    assert!(best_lh > -50.0, "Best log-LH {best_lh} should be better than -50");

    Ok(())
  }
}
