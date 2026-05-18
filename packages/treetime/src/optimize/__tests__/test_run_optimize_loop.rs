#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::run_loop::{ConvergenceReason, run_optimize_loop};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use parking_lot::RwLock;
  use statrs::function::factorial::ln_factorial;
  use std::sync::Arc;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  fn manual_indel_count_on_edge(
    sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
    edge_key: GraphEdgeKey,
  ) -> usize {
    sparse_partitions
      .iter()
      .map(|partition| {
        partition
          .read_arc()
          .edges
          .get(&edge_key)
          .map_or(0, |edge| edge.indels.len())
      })
      .sum()
  }

  fn manual_poisson_indel_log_lh(k: usize, mu: f64, t: f64) -> f64 {
    if k > 0 && t <= 0.0 {
      return f64::NEG_INFINITY;
    }
    if mu == 0.0 {
      return if k == 0 { 0.0 } else { f64::NEG_INFINITY };
    }
    if k == 0 {
      return -mu * t;
    }

    let lambda = mu * t;
    (k as f64) * lambda.ln() - lambda - ln_factorial(k as u64)
  }

  fn manual_total_indel_log_lh(
    graph: &GraphAncestral,
    sparse_partitions: &[Arc<RwLock<PartitionMarginalSparse>>],
  ) -> f64 {
    let total_indels: usize = graph
      .get_edges()
      .iter()
      .map(|edge_ref| manual_indel_count_on_edge(sparse_partitions, edge_ref.read_arc().key()))
      .sum();
    let total_branch_length: f64 = graph
      .get_edges()
      .iter()
      .map(|edge_ref| edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .sum();
    let indel_rate = if total_indels > 0 && total_branch_length > 0.0 {
      total_indels as f64 / total_branch_length
    } else {
      0.0
    };

    graph
      .get_edges()
      .iter()
      .map(|edge_ref| {
        let edge_key = edge_ref.read_arc().key();
        let branch_length = edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0);
        let indel_count = manual_indel_count_on_edge(sparse_partitions, edge_key);
        manual_poisson_indel_log_lh(indel_count, indel_rate, branch_length)
      })
      .sum()
  }

  // Each executed iteration records exactly one log-likelihood entry. With
  // `dp = 0.0` and `damping = 0.75` the convergence and oscillation checks
  // never fire, but the worsened check may stop the loop early once the
  // likelihood begins to decrease.
  #[test]
  fn test_run_optimize_loop_records_lh_history() -> Result<(), Report> {
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
      false,
    )?;

    // One entry per executed iteration, regardless of how the loop stopped.
    assert!(result.lh_history.len() >= 1);
    assert!(result.lh_history.len() <= max_iter);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_records_joint_likelihood_with_sparse_indels() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let first_edge_key = graph.get_edges()[0].read_arc().key();
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.1));
    sparse_partitions[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];

    let sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let dense_lh = update_marginal(&graph, &dense_partitions)?;
    let indel_lh = manual_total_indel_log_lh(&graph, &sparse_partitions);
    let expected_total_lh = sparse_lh + dense_lh + indel_lh;

    let result = run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
    )?;

    assert_eq!(result.lh_history.len(), 1);
    assert_abs_diff_eq!(result.lh_history[0], expected_total_lh, epsilon = 1e-10);
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
      false,
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
      false,
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
      false,
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
    let first_edge_key = graph.get_edges()[0].read_arc().key();
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.1));
    sparse_partitions[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];

    let initial_sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let initial_dense_lh = update_marginal(&graph, &dense_partitions)?;
    let initial_lh = initial_sparse_lh + initial_dense_lh + manual_total_indel_log_lh(&graph, &sparse_partitions);

    run_optimize_loop(
      &mut graph,
      &sparse_partitions,
      &dense_partitions,
      &mixed_partitions,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
    )?;

    let final_sparse_lh = update_marginal(&graph, &sparse_partitions)?;
    let final_dense_lh = update_marginal(&graph, &dense_partitions)?;
    let final_lh = final_sparse_lh + final_dense_lh + manual_total_indel_log_lh(&graph, &sparse_partitions);

    assert!(
      final_lh >= initial_lh,
      "Loop regressed likelihood: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }
}
