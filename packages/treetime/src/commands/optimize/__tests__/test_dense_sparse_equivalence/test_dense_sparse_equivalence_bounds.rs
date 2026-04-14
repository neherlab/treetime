#[cfg(test)]
mod tests {
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::optimize_unified::run_optimize_mixed;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_dense_sparse_equivalence_support::tests::{
    TREE_NEWICK, gap_free_alignment, get_branch_lengths, setup_dense_only, setup_sparse_only,
  };

  #[test]
  fn test_dense_sparse_log_lh_bounded_difference_after_optimization() -> Result<(), Report> {
    let aln = gap_free_alignment()?;

    // Run dense-only optimization
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph_dense, &dense_partitions, BranchOptMethod::Newton)?;
      update_marginal(&graph_dense, &dense_partitions)?;
    }

    let log_lh_dense = update_marginal(&graph_dense, &dense_partitions)?;

    // Run sparse-only optimization
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph_sparse, &sparse_partitions, BranchOptMethod::Newton)?;
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

    // Dense and sparse should converge to similar values.
    // Both use the unified optimizer; differences arise from coefficient
    // representation (per-position dense vs multiplicity-weighted sparse)
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
      run_optimize_mixed(&graph_dense, &dense_partitions, BranchOptMethod::Newton)?;
      update_marginal(&graph_dense, &dense_partitions)?;
    }

    let branch_lengths_dense = get_branch_lengths(&graph_dense);

    // Run sparse-only optimization
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;

    for _ in 0..10 {
      run_optimize_mixed(&graph_sparse, &sparse_partitions, BranchOptMethod::Newton)?;
      update_marginal(&graph_sparse, &sparse_partitions)?;
    }

    let branch_lengths_sparse = get_branch_lengths(&graph_sparse);

    // Both modes should produce same number of edges
    assert_eq!(branch_lengths_dense.len(), branch_lengths_sparse.len());

    // All branch lengths should be valid and bounded
    assert!(
      branch_lengths_dense
        .iter()
        .all(|bl| bl.is_finite() && *bl >= 0.0 && *bl < 10.0),
      "All dense branch lengths should be finite, non-negative, and < 10: {branch_lengths_dense:?}"
    );
    assert!(
      branch_lengths_sparse
        .iter()
        .all(|bl| bl.is_finite() && *bl >= 0.0 && *bl < 10.0),
      "All sparse branch lengths should be finite, non-negative, and < 10: {branch_lengths_sparse:?}"
    );

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
    let max_diff = branch_lengths_dense
      .iter()
      .zip(branch_lengths_sparse.iter())
      .map(|(d, s)| (d - s).abs())
      .fold(0.0_f64, f64::max);
    assert!(
      max_diff < 0.05,
      "Branch lengths should be similar: max_diff={max_diff}, dense={branch_lengths_dense:?}, sparse={branch_lengths_sparse:?}"
    );

    Ok(())
  }
}
