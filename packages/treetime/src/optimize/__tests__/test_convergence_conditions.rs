#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::iteration::{DAMPING_FLOOR, apply_damping};
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{
    ConvergenceReason, marginal_update_dense, marginal_update_sparse, run_optimize_loop,
  };
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use num_traits::pow::pow;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  // At very high iteration counts, the exponential damping factor decays below the floor.
  // The floor ensures the old-value weight never drops below DAMPING_FLOOR.
  #[rustfmt::skip]
  #[rstest]
  #[case::iter_100(  100, DAMPING_FLOOR)]
  #[case::iter_500(  500, DAMPING_FLOOR)]
  #[case::iter_1000(1000, DAMPING_FLOOR)]
  #[trace]
  fn test_convergence_conditions_damping_floor_at_high_iteration(#[case] iteration: usize,
    #[case] expected_old_weight: f64,
  ) -> Result<(), Report> {
    let damping = 0.75;
    let NwkParse { graph, names, branch_lengths, .. } = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;
    let graph: Graph = graph;
    let old_bls = branch_lengths;

    // Set all "optimized" branch lengths to zero.
    let mut bls: BTreeMap<GraphEdgeKey, Option<f64>> = old_bls.keys().map(|&key| (key, Some(0.0))).collect();

    apply_damping(&mut bls, &old_bls, damping, iteration);

    // bl = 0.0 * (1 - old_weight) + 1.0 * old_weight = old_weight
    for bl in bls.values() {
      assert_abs_diff_eq!(bl.unwrap(), expected_old_weight, epsilon = 1e-15);
    }
    Ok(())
  }

  // Verify that below the crossover iteration, the exponential decay is used (not the floor).
  #[test]
  fn test_convergence_conditions_damping_uses_exponential_below_crossover() -> Result<(), Report> {
    let damping = 0.75;
    let iteration = 5; // 0.75^6 = 0.178 >> DAMPING_FLOOR
    let expected_old_weight = pow(damping, iteration + 1);
    assert!(expected_old_weight > DAMPING_FLOOR);

    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;

    let graph: Graph = graph;
    let old_bls = branch_lengths;

    let mut bls: BTreeMap<GraphEdgeKey, Option<f64>> = old_bls.keys().map(|&key| (key, Some(0.0))).collect();

    apply_damping(&mut bls, &old_bls, damping, iteration);

    for bl in bls.values() {
      assert_abs_diff_eq!(bl.unwrap(), expected_old_weight, epsilon = 1e-15);
    }
    Ok(())
  }

  // The convergence check fires when successive likelihoods are within dp.
  // On a toy tree with damping, the loop should converge within a few iterations.
  #[test]
  fn test_convergence_conditions_converged_reason() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut graph: Graph = graph;
    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_6 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      20,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_6,
    )?;

    let (iter, reason) = result.stopped_at.expect("loop should have stopped");
    assert!(
      reason == ConvergenceReason::Converged || reason == ConvergenceReason::Oscillating,
      "Expected Converged or Oscillating on toy tree, got {reason:?} at iteration {iter}"
    );
    Ok(())
  }

  // The worsened condition fires when the likelihood decreases after the peak.
  // On an undamped toy tree, oscillation causes the worsened condition to fire.
  #[test]
  fn test_convergence_conditions_worsened_reverts_to_best() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut graph: Graph = graph;
    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    // Undamped with dp=0 (convergence/oscillation checks never fire) forces the
    // worsened condition to be the only active stopping criterion.
    let names_tt_5 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      50,
      0.0,
      0.0,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_5,
    )?;

    match result.stopped_at {
      Some((iter, ConvergenceReason::Worsened)) => {
        assert!(iter >= 2, "Worsened should not fire before iteration 2, got {iter}");
        // The best LH should be the maximum in the history (the worsened condition
        // restores branch lengths from the best iteration).
        let best_lh = result
          .lh_history
          .iter()
          .map(|log_lh| log_lh.value())
          .fold(f64::NEG_INFINITY, f64::max);
        // The iteration that triggered worsened must have a lower LH than the best.
        let trigger_lh = result.lh_history[iter].value();
        assert!(
          trigger_lh < best_lh,
          "Trigger LH ({trigger_lh:.6}) should be less than best LH ({best_lh:.6})"
        );
      },
      other => panic!("Expected Worsened on undamped toy tree, got {other:?}"),
    }
    Ok(())
  }

  // Rollback validation (T1.4 gate). On a worsening iteration the loop restores the best-seen
  // branch-length map and recomputes the partitions from it, so the returned map must reproduce
  // the best likelihood. With `no_indels = true` the total log-likelihood is exactly the sum of
  // the sparse and dense marginal passes, so recomputing those two from the returned map must
  // match the maximum recorded in `lh_history` (the best likelihood the loop retained).
  // Oracle: the rollback contract in `run_optimize_loop` (restore best map, recompute marginal).
  #[test]
  fn test_convergence_conditions_worsened_rollback_reproduces_best_lh() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut graph: Graph = graph;
    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_4 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      50,
      0.0,
      0.0,
      BranchOptMethod::BrentSqrt,
      true,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_4,
    )?;

    let (_iter, reason) = result.stopped_at.expect("loop should stop");
    assert_eq!(reason, ConvergenceReason::Worsened);

    let best_lh = result
      .lh_history
      .iter()
      .map(|log_lh| log_lh.value())
      .fold(f64::NEG_INFINITY, f64::max);

    // Recompute the marginal likelihood from the returned (rolled-back) branch-length map.
    let marginal_bl = profile_branch_lengths(&result.branch_lengths);
    let sparse_lh = marginal_update_sparse(&graph, &marginal_bl, &mut sparse_partitions)?.value();
    let dense_lh = marginal_update_dense(&graph, &marginal_bl, &mut dense_partitions)?.value();
    assert_abs_diff_eq!(sparse_lh + dense_lh, best_lh, epsilon = 1e-9);
    Ok(())
  }

  // The oscillation detection fires when |LH[i] - LH[i-2]| < dp.
  // Use a moderate dp that catches the 2-cycle amplitude on the toy tree.
  #[test]
  fn test_convergence_conditions_oscillation_detection() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut graph: Graph = graph;
    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    // Use damping to prevent the worsened condition from firing, but set dp
    // large enough that the oscillation check catches the 2-cycle.
    let names_tt_3 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      50,
      1.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_3,
    )?;

    let (iter, reason) = result.stopped_at.expect("loop should have stopped");
    // With dp=1.0, either convergence or oscillation should fire early.
    assert!(
      reason == ConvergenceReason::Converged || reason == ConvergenceReason::Oscillating,
      "Expected early stop with dp=1.0, got {reason:?} at iteration {iter}"
    );
    // Should stop well before max_iter
    assert!(iter < 10, "Expected early stop, got iteration {iter}");
    Ok(())
  }

  // Exhausting max_iter without any stopping condition results in stopped_at = None.
  #[test]
  fn test_convergence_conditions_exhausts_max_iter() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut graph: Graph = graph;
    let (mut dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    // Only 2 iterations with dp=0 and damping. The worsened condition requires
    // i >= 2, so with max_iter=2 (iterations 0 and 1) it cannot fire.
    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      &mut sparse_partitions,
      &mut dense_partitions,
      2,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_2,
    )?;

    assert_eq!(result.lh_history.len(), 2);
    assert!(
      result.stopped_at.is_none(),
      "Should exhaust max_iter=2 without stopping"
    );
    Ok(())
  }

  // Dense-only convergence: verify the fix does not regress dense mode.
  #[test]
  fn test_convergence_conditions_dense_only_converges() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut graph: Graph = graph;

    // Use the setup but only dense partitions (sparse empty)
    let (mut dense_partitions, _sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let mut empty_sparse = vec![];

    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      &mut empty_sparse,
      &mut dense_partitions,
      10,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "Dense-only optimization should converge within 10 iterations"
    );
    Ok(())
  }
}
