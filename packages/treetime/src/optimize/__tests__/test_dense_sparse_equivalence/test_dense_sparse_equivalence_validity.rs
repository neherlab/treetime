#[cfg(test)]
mod tests {
  use super::super::test_dense_sparse_equivalence_support::tests::{
    TREE_NEWICK, gap_free_alignment, setup_dense_only, setup_sparse_only,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts};
  use crate::optimize::params::BranchOptMethod;
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;

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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_dense_only(&graph, &names, &aln, &branch_lengths)?;

    let (mut reconstruction, initial_lh) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let initial_lh = initial_lh.value();
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      let total_length = reconstruction.sequence_length();
      let contributions = gather_edge_contributions(&graph, &reconstruction)?;
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
      let lh;
      (reconstruction, lh) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      let lh = lh.value();
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization");
    }

    let (reconstruction, final_lh) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let final_lh = final_lh.value();

    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization should not significantly decrease likelihood: initial={initial_lh}, final={final_lh}"
    );

    for edge in graph.get_edges() {
      let bl = branch_lengths[&edge.key()].expect("branch length must be set on every edge after optimization");
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_sparse_only(&graph, &names, &aln, &branch_lengths)?;

    let (mut reconstruction, initial_lh) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let initial_lh = initial_lh.value();
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      let total_length = reconstruction.sequence_length();
      let contributions = gather_edge_contributions(&graph, &reconstruction)?;
      let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
      let lh;
      (reconstruction, lh) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
      let lh = lh.value();
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization");
    }

    let (reconstruction, final_lh) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;
    let final_lh = final_lh.value();

    assert!(
      final_lh > -100.0 && final_lh < -10.0,
      "Final log-LH {final_lh} should be in range [-100, -10]"
    );

    assert!(
      final_lh >= initial_lh - 1.0,
      "Optimization should not significantly decrease likelihood: initial={initial_lh}, final={final_lh}"
    );

    for edge in graph.get_edges() {
      let bl = branch_lengths[&edge.key()].expect("branch length must be set on every edge after optimization");
      assert!(bl.is_finite(), "Branch length should be finite");
      assert!(bl >= 0.0, "Branch length should be non-negative");
      assert!(bl < 10.0, "Branch length {bl} should be reasonable (< 10)");
    }

    Ok(())
  }
}
