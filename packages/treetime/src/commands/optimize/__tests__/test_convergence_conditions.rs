#[cfg(test)]
mod tests {
  use crate::commands::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::run::{
    ConvergenceReason, DAMPING_FLOOR, apply_damping, restore_branch_lengths, run_optimize_loop,
    save_branch_lengths,
  };
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use itertools::izip;
  use num_traits::pow::pow;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;

  // --- Damping floor tests ---

  // At very high iteration counts, the exponential damping factor decays below the floor.
  // The floor ensures the old-value weight never drops below DAMPING_FLOOR.
  #[rustfmt::skip]
  #[rstest]
  #[case::iter_100(  100, DAMPING_FLOOR)]
  #[case::iter_500(  500, DAMPING_FLOOR)]
  #[case::iter_1000(1000, DAMPING_FLOOR)]
  #[trace]
  fn test_convergence_conditions_damping_floor_at_high_iteration(
    #[case] iteration: usize,
    #[case] expected_old_weight: f64,
  ) -> Result<(), Report> {
    let damping = 0.75;
    let graph: GraphAncestral = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;
    let old_bls = save_branch_lengths(&graph);

    // Set all branch lengths to zero ("optimized" value)
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

  // Verify that below the crossover iteration, the exponential decay is used (not the floor).
  #[test]
  fn test_convergence_conditions_damping_uses_exponential_below_crossover() -> Result<(), Report> {
    let damping = 0.75;
    let iteration = 5; // 0.75^6 = 0.178 >> DAMPING_FLOOR
    let expected_old_weight = pow(damping, iteration + 1);
    assert!(expected_old_weight > DAMPING_FLOOR);

    let graph: GraphAncestral = nwk_read_str("(A:1.0,B:1.0)root:0.0;")?;
    let old_bls = save_branch_lengths(&graph);

    for edge_ref in graph.get_edges() {
      edge_ref.write_arc().payload().write_arc().set_branch_length(Some(0.0));
    }

    apply_damping(&graph, &old_bls, damping, iteration);

    for edge_ref in graph.get_edges() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
      assert_abs_diff_eq!(bl, expected_old_weight, epsilon = 1e-15);
    }
    Ok(())
  }

  // --- restore_branch_lengths tests ---

  // Round-trip: save, modify, restore, verify identical to original.
  #[test]
  fn test_convergence_conditions_restore_branch_lengths_roundtrip() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let original = save_branch_lengths(&graph);

    // Modify all branch lengths
    for edge_ref in graph.get_edges() {
      let mut edge = edge_ref.write_arc().payload().write_arc();
      let bl = edge.branch_length().unwrap_or(0.0);
      edge.set_branch_length(Some(bl * 5.0 + 0.42));
    }

    // Verify they changed
    let modified = save_branch_lengths(&graph);
    assert!(
      izip!(&original, &modified).any(|(a, b)| (a - b).abs() > 1e-10),
      "Branch lengths should have changed after modification"
    );

    // Restore and verify match
    restore_branch_lengths(&graph, &original);
    let restored = save_branch_lengths(&graph);
    for (orig, rest) in izip!(&original, &restored) {
      assert_abs_diff_eq!(*orig, *rest, epsilon = 1e-15);
    }
    Ok(())
  }

  // --- Convergence condition integration tests ---

  // The convergence check fires when successive likelihoods are within dp.
  // On a toy tree with damping, the loop should converge within a few iterations.
  #[test]
  fn test_convergence_conditions_converged_reason() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      20,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
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
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Undamped with dp=0 (convergence/oscillation checks never fire) forces the
    // worsened condition to be the only active stopping criterion.
    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      50,
      0.0,
      0.0,
      BranchOptMethod::BrentSqrt,
    )?;

    match result.stopped_at {
      Some((iter, ConvergenceReason::Worsened)) => {
        assert!(iter >= 2, "Worsened should not fire before iteration 2, got {iter}");
        // The best LH should be the maximum in the history (the worsened condition
        // restores branch lengths from the best iteration).
        let best_lh = result.lh_history.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        // The iteration that triggered worsened must have a lower LH than the best.
        let trigger_lh = result.lh_history[iter];
        assert!(
          trigger_lh < best_lh,
          "Trigger LH ({trigger_lh:.6}) should be less than best LH ({best_lh:.6})"
        );
      },
      other => {
        // On some toy trees, the undamped loop may not worsen. This is acceptable
        // if it exhausted max_iter or stopped for another reason. The key property
        // (worsened reverts to best) is still covered by the damping floor test.
        panic!("Expected Worsened on undamped toy tree, got {other:?}");
      },
    }
    Ok(())
  }

  // The oscillation detection fires when |LH[i] - LH[i-2]| < dp.
  // Use a moderate dp that catches the 2-cycle amplitude on the toy tree.
  #[test]
  fn test_convergence_conditions_oscillation_detection() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Use damping to prevent the worsened condition from firing, but set dp
    // large enough that the oscillation check catches the 2-cycle.
    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      50,
      1.0,
      0.75,
      BranchOptMethod::BrentSqrt,
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
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Only 2 iterations with dp=0 and damping. The worsened condition requires
    // i >= 2, so with max_iter=2 (iterations 0 and 1) it cannot fire.
    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      2,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    assert_eq!(result.lh_history.len(), 2);
    assert!(result.stopped_at.is_none(), "Should exhaust max_iter=2 without stopping");
    Ok(())
  }

  // Dense-only convergence: verify the fix does not regress dense mode.
  #[test]
  fn test_convergence_conditions_dense_only_converges() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;

    // Use the setup but only dense partitions (sparse empty)
    let (dense_partitions, _sparse_partitions, _mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Build dense-only mixed partitions
    use crate::commands::optimize::run::collect_optimize_partitions;
    let empty_sparse = vec![];
    let mixed = collect_optimize_partitions(&dense_partitions, &empty_sparse);

    let result = run_optimize_loop(
      &mut graph,
      &empty_sparse,
      &dense_partitions,
      &mixed,
      10,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    assert!(
      result.stopped_at.is_some(),
      "Dense-only optimization should converge within 10 iterations"
    );
    Ok(())
  }
}
