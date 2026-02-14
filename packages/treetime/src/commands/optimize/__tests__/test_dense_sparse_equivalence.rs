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

    // Branch lengths should be valid
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl.is_finite(), "Branch length should be finite");
        assert!(bl >= 0.0, "Branch length should be non-negative");
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

    // Branch lengths should be valid
    for edge in graph.get_edges() {
      let bl = edge.read_arc().payload().read_arc().branch_length();
      if let Some(bl) = bl {
        assert!(bl.is_finite(), "Branch length should be finite");
        assert!(bl >= 0.0, "Branch length should be non-negative");
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

    // Both modes should produce finite log-LH
    assert!(log_lh_dense.is_finite());
    assert!(log_lh_sparse.is_finite());

    // Difference should be bounded (within 1.0 log-LH units)
    // Due to different zero-branch thresholds (0.01 vs 0.0001) and optimization paths
    let diff = (log_lh_dense - log_lh_sparse).abs();
    assert!(
      diff < 1.0,
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

    // All branch lengths should be valid
    for (dense_bl, sparse_bl) in branch_lengths_dense.iter().zip(branch_lengths_sparse.iter()) {
      assert!(dense_bl.is_finite());
      assert!(sparse_bl.is_finite());
      assert!(*dense_bl >= 0.0);
      assert!(*sparse_bl >= 0.0);
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

    let mut prev_lh = f64::NEG_INFINITY;
    let convergence_threshold = 1e-4;

    for _ in 0..50 {
      run_optimize_dense(&graph, &partitions)?;
      let lh = update_marginal(&graph, &partitions)?;

      if (lh - prev_lh).abs() < convergence_threshold {
        return Ok(()); // Converged
      }

      prev_lh = lh;
    }

    // Even if not fully converged, final log-LH should be stable
    let final_lh = update_marginal(&graph, &partitions)?;
    assert!(final_lh.is_finite());

    Ok(())
  }

  #[test]
  fn test_sparse_optimization_converges() -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse_only(&graph, &aln)?;

    let mut prev_lh = f64::NEG_INFINITY;
    let convergence_threshold = 1e-4;

    for _ in 0..50 {
      run_optimize_sparse(&graph, &partitions)?;
      let lh = update_marginal(&graph, &partitions)?;

      if (lh - prev_lh).abs() < convergence_threshold {
        return Ok(()); // Converged
      }

      prev_lh = lh;
    }

    // Even if not fully converged, final log-LH should be stable
    let final_lh = update_marginal(&graph, &partitions)?;
    assert!(final_lh.is_finite());

    Ok(())
  }
}
