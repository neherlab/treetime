#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let partitions = setup_dense_only(&graph, &names, &aln, &branch_lengths)?;

    let (mut partitions, initial_lh) = marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
    let initial_lh = initial_lh.value();
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      let total_length = total_sequence_length(&partitions, &[]);
      let contributions = gather_edge_contributions(&graph, &partitions, &[])?;
      let indel_counts = gather_edge_indel_counts(&graph, &partitions, &[]);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
      let lh;
      (partitions, lh) = marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
      let lh = lh.value();
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization");
    }

    let (partitions, final_lh) = marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
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
    let graph: Graph = graph;
    let partitions = setup_sparse_only(&graph, &names, &aln, &branch_lengths)?;

    let (mut partitions, initial_lh) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
    let initial_lh = initial_lh.value();
    assert!(initial_lh.is_finite(), "Initial log-LH should be finite");

    for _ in 0..10 {
      let total_length = total_sequence_length(&[], &partitions);
      let contributions = gather_edge_contributions(&graph, &[], &partitions)?;
      let indel_counts = gather_edge_indel_counts(&graph, &[], &partitions);
      run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;
      let lh;
      (partitions, lh) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
      let lh = lh.value();
      assert!(lh.is_finite(), "Log-LH should remain finite during optimization");
    }

    let (partitions, final_lh) = marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), partitions)?;
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
