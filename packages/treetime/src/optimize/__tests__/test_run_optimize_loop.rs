#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::params::{BranchOptMethod, TopologyOps};
  use crate::optimize::run_loop::{ConvergenceReason, run_optimize_loop};
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use statrs::function::factorial::ln_factorial;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  fn manual_indel_count_on_edge(sparse_partitions: &[SparseReconstruction], edge_key: GraphEdgeKey) -> usize {
    sparse_partitions
      .iter()
      .map(|partition| {
        partition
          .partition
          .obs_edges
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
    graph: &Graph,
    sparse_partitions: &[SparseReconstruction],
    branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> f64 {
    let total_indels: usize = graph
      .get_edges()
      .map(|edge_ref| manual_indel_count_on_edge(sparse_partitions, edge_ref.key()))
      .sum();
    let total_branch_length: f64 = graph
      .get_edges()
      .map(|edge_ref| branch_lengths.get(&edge_ref.key()).copied().flatten().unwrap_or(0.0))
      .sum();
    let indel_rate = if total_indels > 0 && total_branch_length > 0.0 {
      total_indels as f64 / total_branch_length
    } else {
      0.0
    };

    graph
      .get_edges()
      .map(|edge_ref| {
        let edge_key = edge_ref.key();
        let branch_length = branch_lengths.get(&edge_key).copied().flatten().unwrap_or(0.0);
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 5;
    let names_tt_6 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      max_iter,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_6,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    // One entry per executed iteration, regardless of how the loop stopped.
    assert!(result.lh_history.len() >= 1);
    assert!(result.lh_history.len() <= max_iter);
    Ok(())
  }

  #[test]
  fn test_run_optimize_loop_records_joint_likelihood_with_sparse_indels() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.1));
    sparse_partitions[0]
      .partition
      .obs_edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let (sparse_partitions, sparse_lh) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
    let sparse_lh = sparse_lh.value();
    let (dense_partitions, dense_lh) =
      marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;
    let dense_lh = dense_lh.value();
    let indel_lh = manual_total_indel_log_lh(&graph, &sparse_partitions, &branch_lengths);
    let expected_total_lh = sparse_lh + dense_lh + indel_lh;

    let names_tt_5 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_5,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert_eq!(result.lh_history.len(), 1);
    assert_abs_diff_eq!(result.lh_history[0].value(), expected_total_lh, epsilon = 1e-10);
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let max_iter = 50;
    let dp = f64::INFINITY;

    let names_tt_4 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      max_iter,
      dp,
      0.0,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_4,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    assert_eq!(result.stopped_at, Some((1, ConvergenceReason::Converged)));
    assert_eq!(result.lh_history.len(), 2);
    Ok(())
  }

  // With `max_iter = 0`, the loop body never runs.
  #[test]
  fn test_run_optimize_loop_zero_max_iter_is_noop() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_3 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      0,
      1e-2,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_3,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let names_tt_2 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_2,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;

    for (i, lh) in result.lh_history.iter().map(|log_lh| log_lh.value()).enumerate() {
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
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let mut graph: Graph = graph;
    let (dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.1));
    sparse_partitions[0]
      .partition
      .obs_edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let (sparse_partitions, initial_sparse_lh) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
    let initial_sparse_lh = initial_sparse_lh.value();
    let (dense_partitions, initial_dense_lh) =
      marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;
    let initial_dense_lh = initial_dense_lh.value();
    let initial_lh =
      initial_sparse_lh + initial_dense_lh + manual_total_indel_log_lh(&graph, &sparse_partitions, &branch_lengths);

    let names_tt_1 = names.clone();
    let result = run_optimize_loop(
      &mut graph,
      sparse_partitions,
      dense_partitions,
      10,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths,
      &names_tt_1,
    )?;
    let sparse_partitions = result.sparse_partitions;
    let dense_partitions = result.dense_partitions;
    let branch_lengths = result.branch_lengths;

    let (sparse_partitions, final_sparse_lh) =
      marginal_update_sparse(&graph, &branch_lengths_or_zero(&branch_lengths), sparse_partitions)?;
    let final_sparse_lh = final_sparse_lh.value();
    let (dense_partitions, final_dense_lh) =
      marginal_update_dense(&graph, &branch_lengths_or_zero(&branch_lengths), dense_partitions)?;
    let final_dense_lh = final_dense_lh.value();
    let final_lh =
      final_sparse_lh + final_dense_lh + manual_total_indel_log_lh(&graph, &sparse_partitions, &branch_lengths);

    assert!(
      final_lh >= initial_lh,
      "Loop regressed likelihood: {initial_lh:.6} -> {final_lh:.6}"
    );
    Ok(())
  }
}
