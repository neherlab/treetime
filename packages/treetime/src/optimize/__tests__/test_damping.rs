#[cfg(test)]
mod tests {
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, compute_total_lh, setup_partitions, simple_alignment,
  };
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::iteration::{apply_damping, save_branch_lengths};
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::partition::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::izip;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_save_branch_lengths_captures_all_edges() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let saved = save_branch_lengths(&graph);
    let edges = graph.get_edges();
    assert_eq!(saved.len(), edges.len());
    for (bl, edge_ref) in izip!(&saved, &edges) {
      let expected = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
      assert_abs_diff_eq!(*bl, expected, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_zero_is_noop() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let original = save_branch_lengths(&graph);

    // Manually change branch lengths to simulate optimization
    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();
      let bl = edge.branch_length().unwrap_or(0.0);
      edge.set_branch_length(Some(bl * 2.0));
    }

    let after_optim = save_branch_lengths(&graph);
    apply_damping(&graph, &original, 0.0, 0);
    let after_damping = save_branch_lengths(&graph);

    // damping=0.0 should leave the optimized values untouched
    for (damped, optimized) in izip!(&after_damping, &after_optim) {
      assert_abs_diff_eq!(*damped, *optimized, epsilon = 1e-15);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::iter_0( 0, 0.750)]
  #[case::iter_1( 1, 0.5625)]
  #[case::iter_2( 2, 0.421875)]
  #[case::iter_4( 4, 0.2373046875)]
  #[case::iter_9( 9, 0.056313514709472656)]
  #[trace]
  fn test_apply_damping_weights_match_v0(#[case] iteration: usize, #[case] expected_old_weight: f64) -> Result<(), Report> {
    let damping = 0.75;
    let graph: GraphAncestral = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;
    let old_bls = save_branch_lengths(&graph);

    // Set all branch lengths to a known "optimized" value
    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }

    apply_damping(&graph, &old_bls, damping, iteration);

    // bl = 0.0 * (1 - old_weight) + 1.0 * old_weight = old_weight
    for edge_ref in graph.get_edges() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
      assert_abs_diff_eq!(bl, expected_old_weight, epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_apply_damping_blends_correctly() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.2)root:0.0;")?;
    let old_bls = save_branch_lengths(&graph);

    // Set "optimized" branch lengths
    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();
      let bl = edge.branch_length().unwrap_or(0.0);
      edge.set_branch_length(Some(bl * 3.0));
    }

    apply_damping(&graph, &old_bls, 0.75, 0);

    // At iteration 0, damping_factor = 0.75, new_weight = 0.25
    // Edge A: 0.3 * 0.25 + 0.1 * 0.75 = 0.075 + 0.075 = 0.15
    // Edge B: 0.6 * 0.25 + 0.2 * 0.75 = 0.15 + 0.15 = 0.30
    // Epsilon accounts for Newick float parsing roundtrip
    let damped = save_branch_lengths(&graph);
    assert_abs_diff_eq!(damped[0], 0.15, epsilon = 1e-8);
    assert_abs_diff_eq!(damped[1], 0.30, epsilon = 1e-8);
    Ok(())
  }

  #[test]
  fn test_apply_damping_new_weight_increases_with_iteration() -> Result<(), Report> {
    let damping = 0.75;
    let old_bl = 1.0;
    let optimized_bl = 2.0;

    let mut prev_damped = old_bl;
    for iteration in 0..10 {
      let graph: GraphAncestral = nwk_read_str("(A:1.0)root:0.0;")?;
      let old_bls = save_branch_lengths(&graph);
      graph
        .get_edges()
        .first()
        .unwrap()
        .write_arc()
        .payload()
        .write_arc()
        .set_branch_length(Some(optimized_bl));

      apply_damping(&graph, &old_bls, damping, iteration);

      let damped = graph
        .get_edges()
        .first()
        .unwrap()
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap_or(0.0);

      // Each subsequent iteration should give more weight to the optimized value
      assert!(
        damped > prev_damped || iteration == 0,
        "Iteration {iteration}: damped {damped} should be > previous {prev_damped}"
      );
      prev_damped = damped;
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
  fn test_damped_optimization_converges(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let max_iter = 10;
    let damping = 0.75;
    let dp = 0.1;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      max_iter,
      dp,
      damping,
      method,
      false,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "Damped optimization did not stop within {max_iter} iterations"
    );

    // Final log-likelihood must be within a tight range around the observed fixed point.
    // The toy tree (4 leaves, 16 sites, JC69) converges near -72.41. Use the post-loop
    // marginal pass so `final_lh` reflects the state after the last branch-length update,
    // not the pre-update measurement recorded in `lh_history`.
    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(
      final_lh > -73.0 && final_lh < -72.0,
      "Final log-lh {final_lh:.6} outside expected range (-73.0, -72.0)"
    );

    // The three-condition convergence check (converged, oscillating, worsened) detects
    // the 2-cycle on this toy tree. The loop stops via the oscillating or converged
    // condition before sign flips accumulate.

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
  fn test_damped_optimization_does_not_regress(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let initial_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;

    // Force all 10 iterations (never break on convergence) so the non-regression check
    // exercises the full damping trajectory rather than possibly stopping after two
    // near-identical likelihoods.
    let dp = 0.0;
    run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      10,
      dp,
      0.75,
      method,
      false,
    )?;

    // Strict non-regression: damped optimization must not degrade likelihood.
    // Damping blends new and old branch lengths as a convex combination,
    // so overall likelihood should improve or hold steady.
    let final_lh = compute_total_lh(&graph, &dense_partitions, &sparse_partitions)?;
    assert!(
      final_lh >= initial_lh,
      "Damped optimization regressed: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }
}
