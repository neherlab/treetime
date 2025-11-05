use crate::commands::optimize::optimize_unified::{OptimizationContribution, evaluate_mixed};
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use eyre::Report;
use ndarray::{Array1, s};
use std::sync::Arc;

pub fn collect_edge_contributions(
  edge_key: GraphEdgeKey,
  dense_partitions: &[&PartitionMarginalDense],
  sparse_partitions: &[&PartitionMarginalSparse],
) -> Result<Vec<OptimizationContribution>, Report> {
  let mut contributions = Vec::with_capacity(dense_partitions.len() + sparse_partitions.len());

  for partition in dense_partitions {
    contributions.push(OptimizationContribution::from_dense(edge_key, partition));
  }

  for partition in sparse_partitions {
    contributions.push(OptimizationContribution::from_sparse(edge_key, partition)?);
  }

  Ok(contributions)
}

pub fn compute_branch_length_distribution(
  contributions: &[OptimizationContribution],
  current_branch_length: f64,
  one_mutation: f64,
  n_grid_points: usize,
  clock_rate: f64,
) -> Result<Arc<Distribution>, Report> {
  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points);

  let log_prob = grid.mapv(|branch_len| {
    let metrics = evaluate_mixed(contributions, branch_len);
    metrics.log_lh
  });

  let max_log_lh = log_prob.fold(f64::NEG_INFINITY, |a, &b| a.max(b));

  let normalized_prob = log_prob.mapv(|lh| (lh - max_log_lh).exp());

  let time_grid = grid.mapv(|bl| bl / clock_rate);

  if clock_rate < 0.0 {
    let time_grid_rev = time_grid.slice(s![..;-1]).to_owned();
    let normalized_prob_rev = normalized_prob.slice(s![..;-1]).to_owned();
    Distribution::function(time_grid_rev, normalized_prob_rev)
  } else {
    Distribution::function(time_grid, normalized_prob)
  }
  .map(Arc::new)
}

fn create_simple_grid(center: f64, one_mutation: f64, n_points: usize) -> Array1<f64> {
  let min_bl = one_mutation * 0.1;
  let max_bl = f64::max(center * 3.0, one_mutation * 10.0);

  Array1::linspace(min_bl, max_bl, n_points)
}
