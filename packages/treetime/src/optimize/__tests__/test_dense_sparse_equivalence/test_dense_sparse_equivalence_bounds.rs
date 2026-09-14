#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_dense_sparse_equivalence_support::tests::{
    TREE_NEWICK, gap_free_alignment, get_branch_lengths, setup_dense_only, setup_sparse_only,
  };

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_dense_sparse_log_lh_bounded_difference_after_optimization(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = gap_free_alignment()?;

    // Run dense-only optimization
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let mut bl_dense = nwk_parsed.branch_lengths;
    let mut dense_partitions = setup_dense_only(&graph_dense, &graph_dense_names, &aln, &bl_dense)?;

    for _ in 0..10 {
      let total_length = total_sequence_length(&dense_partitions, &[]);
      let contributions = gather_edge_contributions(&graph_dense, &dense_partitions, &[])?;
      let indel_counts = gather_edge_indel_counts(&graph_dense, &dense_partitions, &[]);
      run_optimize_mixed(&graph_dense, total_length, &contributions, &indel_counts, method, &mut bl_dense)?;
      (dense_partitions, _) = marginal_update_dense(&graph_dense, &branch_lengths_or_zero(&bl_dense), dense_partitions)?;
    }

    let (dense_partitions, log_lh_dense) = marginal_update_dense(&graph_dense, &branch_lengths_or_zero(&bl_dense), dense_partitions)?;
    let log_lh_dense = log_lh_dense.value();

    // Run sparse-only optimization
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let mut bl_sparse = nwk_parsed.branch_lengths;
    let mut sparse_partitions = setup_sparse_only(&graph_sparse, &graph_sparse_names, &aln, &bl_sparse)?;

    for _ in 0..10 {
      let total_length = total_sequence_length(&[], &sparse_partitions);
      let contributions = gather_edge_contributions(&graph_sparse, &[], &sparse_partitions)?;
      let indel_counts = gather_edge_indel_counts(&graph_sparse, &[], &sparse_partitions);
      run_optimize_mixed(&graph_sparse, total_length, &contributions, &indel_counts, method, &mut bl_sparse)?;
      (sparse_partitions, _) = marginal_update_sparse(&graph_sparse, &branch_lengths_or_zero(&bl_sparse), sparse_partitions)?;
    }

    let (sparse_partitions, log_lh_sparse) = marginal_update_sparse(&graph_sparse, &branch_lengths_or_zero(&bl_sparse), sparse_partitions)?;
    let log_lh_sparse = log_lh_sparse.value();

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

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_dense_sparse_branch_lengths_bounded_difference(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = gap_free_alignment()?;

    // Run dense-only optimization
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let mut bl_dense = nwk_parsed.branch_lengths;
    let mut dense_partitions = setup_dense_only(&graph_dense, &graph_dense_names, &aln, &bl_dense)?;

    for _ in 0..10 {
      let total_length = total_sequence_length(&dense_partitions, &[]);
      let contributions = gather_edge_contributions(&graph_dense, &dense_partitions, &[])?;
      let indel_counts = gather_edge_indel_counts(&graph_dense, &dense_partitions, &[]);
      run_optimize_mixed(&graph_dense, total_length, &contributions, &indel_counts, method, &mut bl_dense)?;
      (dense_partitions, _) = marginal_update_dense(&graph_dense, &branch_lengths_or_zero(&bl_dense), dense_partitions)?;
    }

    let branch_lengths_dense = get_branch_lengths(&graph_dense, &bl_dense);

    // Run sparse-only optimization
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let mut bl_sparse = nwk_parsed.branch_lengths;
    let mut sparse_partitions = setup_sparse_only(&graph_sparse, &graph_sparse_names, &aln, &bl_sparse)?;

    for _ in 0..10 {
      let total_length = total_sequence_length(&[], &sparse_partitions);
      let contributions = gather_edge_contributions(&graph_sparse, &[], &sparse_partitions)?;
      let indel_counts = gather_edge_indel_counts(&graph_sparse, &[], &sparse_partitions);
      run_optimize_mixed(&graph_sparse, total_length, &contributions, &indel_counts, method, &mut bl_sparse)?;
      (sparse_partitions, _) = marginal_update_sparse(&graph_sparse, &branch_lengths_or_zero(&bl_sparse), sparse_partitions)?;
    }

    let branch_lengths_sparse = get_branch_lengths(&graph_sparse, &bl_sparse);

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
