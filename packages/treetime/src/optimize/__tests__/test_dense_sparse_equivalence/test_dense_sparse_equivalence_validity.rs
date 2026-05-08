#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::optimize::args::BranchOptMethod;
  use crate::optimize::optimize_unified::run_optimize_mixed;
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_dense_sparse_equivalence_support::tests::{
    TREE_NEWICK, gap_free_alignment, setup_dense_only, setup_sparse_only,
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
  fn test_dense_optimization_produces_valid_results(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_dense_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      run_optimize_mixed(&graph, &partitions, method)?;
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
      let bl = edge
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .expect("branch length must be set on every edge after optimization");
      assert!(bl.is_finite(), "Branch length should be finite");
      assert!(bl >= 0.0, "Branch length should be non-negative");
      assert!(bl < 10.0, "Branch length {bl} should be reasonable (< 10)");
    }

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
  fn test_sparse_optimization_produces_valid_results(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = gap_free_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let partitions = setup_sparse_only(&graph, &aln)?;

    let initial_lh = update_marginal(&graph, &partitions)?;
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      run_optimize_mixed(&graph, &partitions, method)?;
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
      let bl = edge
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .expect("branch length must be set on every edge after optimization");
      assert!(bl.is_finite(), "Branch length should be finite");
      assert!(bl >= 0.0, "Branch length should be non-negative");
      assert!(bl < 10.0, "Branch length {bl} should be reasonable (< 10)");
    }

    Ok(())
  }
}
