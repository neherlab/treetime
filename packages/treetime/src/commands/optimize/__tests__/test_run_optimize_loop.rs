#[cfg(test)]
mod tests {
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::commands::optimize::args::BranchOptMethod;
  use crate::commands::optimize::run::{ConvergenceReason, run_optimize_loop};
  use crate::representation::payload::ancestral::GraphAncestral;
  use eyre::Report;
  use treetime_io::nwk::nwk_read_str;

  // Each executed iteration records exactly one log-likelihood entry. With
  // `dp = 0.0` and `damping = 0.75` the convergence and oscillation checks
  // never fire, but the worsened check may stop the loop early once the
  // likelihood begins to decrease.
  #[test]
  fn test_run_optimize_loop_records_lh_per_iteration() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let max_iter = 5;
    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      max_iter,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    // One entry per executed iteration, regardless of how the loop stopped.
    assert!(result.lh_history.len() >= 1);
    assert!(result.lh_history.len() <= max_iter);
    Ok(())
  }

  // When `|ΔLH| < |dp|` the loop must break and report the iteration at which
  // that happened via `stopped_at`.
  //
  // With `dp = INFINITY`, the convergence check fires on iteration 1 (the first
  // iteration with a finite `lh_prev`). Iteration 0 computes
  // `|total_lh - NEG_INFINITY| = INFINITY` and `INFINITY < INFINITY` is false,
  // so it does not trigger.
  #[test]
  fn test_run_optimize_loop_breaks_on_convergence() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let max_iter = 50;
    let dp = f64::INFINITY;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      max_iter,
      dp,
      0.0,
      BranchOptMethod::BrentSqrt,
    )?;

    assert_eq!(result.stopped_at, Some((1, ConvergenceReason::Converged)));
    assert_eq!(result.lh_history.len(), 2);
    Ok(())
  }

  // With `max_iter = 0`, the loop body never runs.
  #[test]
  fn test_run_optimize_loop_zero_max_iter_is_noop() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      0,
      1e-2,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    assert!(result.lh_history.is_empty());
    assert!(result.stopped_at.is_none());
    Ok(())
  }

  // All log-likelihood values recorded during the optimization loop must be finite.
  // This guards against NaN/inf from forward-pass division by zero or degenerate
  // normalization. The defect was: unguarded numerator/divisor in the forward pass
  // produced NaN when the divisor contained zeros, poisoning all subsequent iterations.
  #[test]
  fn test_run_optimize_loop_all_likelihoods_finite() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    for (i, &lh) in result.lh_history.iter().enumerate() {
      assert!(lh.is_finite(), "Iteration {i}: log-likelihood must be finite, got {lh}");
    }
    Ok(())
  }

  // Optimization should improve (or maintain) likelihood relative to the initial
  // state. With damping, the loop converges smoothly. The worsened condition may
  // stop the loop before max_iter, but the final state is always the best observed.
  #[test]
  fn test_run_optimize_loop_improves_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let initial_sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let initial_dense_lh = update_marginal(&graph, &dense_partitions)?;
    let initial_lh = initial_sparse_lh + initial_dense_lh;

    run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
    )?;

    let final_sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let final_dense_lh = update_marginal(&graph, &dense_partitions)?;
    let final_lh = final_sparse_lh + final_dense_lh;

    assert!(
      final_lh >= initial_lh,
      "Loop regressed likelihood: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }
}
