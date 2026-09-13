#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::{OptimizeReadouts, marginal_update_dense, marginal_update_sparse};
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  use super::super::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_partitions, simple_alignment,
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
  fn test_optimization_converges_within_iterations(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let max_iter = 50;

    // Reference: run BrentSqrt (the v0-matching default) on a fresh graph
    // and capture its final undamped log-likelihood. Tying the per-method
    // assertion to the reference ties this test to cross-method agreement
    // rather than to a hand-chosen LH range that would silently drift.
    let lh_ref = {
      let NwkParse { graph: graph_ref, names: graph_ref_names, branch_lengths: mut branch_lengths_ref, .. } = nwk_read_str(TREE_NEWICK)?;
      let (dp_ref, sp_ref) = setup_partitions(&graph_ref, &graph_ref_names, &aln, &mut branch_lengths_ref)?;
      let ro_ref = OptimizeReadouts::new(&dp_ref, &sp_ref);
      let mp_ref = ro_ref.view();
      for _ in 0..max_iter {
        run_optimize_mixed(&graph_ref, &mp_ref, BranchOptMethod::BrentSqrt, &mut branch_lengths_ref)?;
      }
      let (_, _, lh) = compute_total_lh(&graph_ref, dp_ref, sp_ref, &branch_lengths_ref)?;
      lh
    };

    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(TREE_NEWICK)?;

    let graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    for _ in 0..max_iter {
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
    }
    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;

    // The undamped alternating optimization produces a stable 2-cycle. Each
    // method must converge to within 1e-2 of the BrentSqrt reference's final
    // log-likelihood. The 1e-2 tolerance accommodates the 2-cycle
    // (sub-iteration variation) for all six methods on this toy alignment.
    let lh_diff = (final_lh - lh_ref).abs();
    assert!(
      lh_diff < 1e-2,
      "{method:?} final lh {final_lh:.6} differs from BrentSqrt reference {lh_ref:.6} by {lh_diff:.6} > 1e-2"
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
  fn test_optimization_improves_or_maintains_likelihood(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(TREE_NEWICK)?;
    let graph: Graph = graph;

    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let (dense_partitions, sparse_partitions, initial_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;
    // Initial log-lh should be negative (log of probability < 1)
    assert!(initial_lh < 0.0, "Initial log-LH should be negative: {initial_lh}");

    // Run several optimization steps
    for _ in 0..10 {
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
    }

    let (dense_partitions, sparse_partitions, final_lh) = compute_total_lh(&graph, dense_partitions, sparse_partitions, &branch_lengths)?;

    // Strict non-regression: optimization must not degrade likelihood.
    // This test runs pure branch length optimization without marginal reconstruction
    // alternation, so there is no 2-cycle and likelihood should improve monotonically.
    assert!(
      final_lh >= initial_lh,
      "Optimization regressed: {initial_lh:.6} -> {final_lh:.6}"
    );

    // Final log-lh should be in reasonable range for JC69 on this alignment
    assert!(final_lh < 0.0, "Log-LH should be negative: {final_lh}");
    assert!(final_lh > -100.0, "Log-LH unreasonably low: {final_lh}");

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
  fn test_optimization_produces_valid_branch_lengths(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse { graph, names, mut branch_lengths, .. } = nwk_read_str(TREE_NEWICK)?;
    let graph: Graph = graph;

    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    // Collect initial branch lengths
    let initial_total: f64 = graph
      .get_edges()
      .iter()
      .filter_map(|e| branch_lengths[&e.read_arc().key()])
      .sum();

    // Run several optimization iterations
    for _ in 0..10 {
      run_optimize_mixed(&graph, &OptimizeReadouts::new(&dense_partitions, &sparse_partitions).view(), method, &mut branch_lengths)?;
      (dense_partitions, _) = marginal_update_dense(&graph, &profile_branch_lengths(&branch_lengths), dense_partitions)?;
      (sparse_partitions, _) = marginal_update_sparse(&graph, &profile_branch_lengths(&branch_lengths), sparse_partitions)?;
    }

    // Verify all branch lengths are in valid range
    let mut final_total = 0.0;
    for edge in graph.get_edges() {
      let edge = edge.read_arc();
      let branch_length = branch_lengths[&edge.key()];
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
}
