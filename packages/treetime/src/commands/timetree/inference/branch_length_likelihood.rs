use crate::commands::optimize::optimize_unified::{OptimizationContribution, evaluate_mixed_log_lh_only};
use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;

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
  debug_assert!(clock_rate >= 0.0);

  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points);

  let log_lh = grid.mapv(|branch_len| evaluate_mixed_log_lh_only(contributions, branch_len));
  let max_log_lh = log_lh.max()?;

  let normalized_prob = (&log_lh - *max_log_lh).exp();

  // Convert branch length grid to time grid by dividing by clock rate
  let time_min = grid[0] / clock_rate;
  let time_max = grid[grid.len() - 1] / clock_rate;

  let distribution_fn = DistributionFunction::from_range_values((time_min, time_max), normalized_prob)?;
  Ok(Arc::new(Distribution::Function(distribution_fn)))
}

fn create_simple_grid(center: f64, one_mutation: f64, n_points: usize) -> Array1<f64> {
  let min_bl = one_mutation * 0.1;
  let max_bl = f64::max(center * 3.0, one_mutation * 10.0);
  Array1::linspace(min_bl, max_bl, n_points)
}
