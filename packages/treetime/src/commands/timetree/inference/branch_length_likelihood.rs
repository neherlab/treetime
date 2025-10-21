use crate::commands::optimize::optimize_unified::{OptimizationContribution, evaluate_mixed};
use crate::distribution::distribution::Distribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use eyre::Report;
use ndarray::Array1;
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
) -> Result<Arc<Distribution>, Report> {
  let grid = create_simple_grid(current_branch_length, one_mutation, n_grid_points);

  let log_prob: Vec<f64> = grid
    .iter()
    .map(|&bl| {
      let metrics = evaluate_mixed(contributions, bl);
      metrics.log_lh
    })
    .collect();

  let max_log_lh = log_prob.iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let normalized_log_prob: Vec<f64> = log_prob.iter().map(|&lh| lh - max_log_lh).collect();

  Distribution::function(Array1::from_vec(grid), Array1::from_vec(normalized_log_prob)).map(Arc::new)
}

fn create_simple_grid(center: f64, one_mutation: f64, n_points: usize) -> Vec<f64> {
  let min_bl = one_mutation * 0.1;
  let max_bl = f64::max(center * 3.0, one_mutation * 10.0);

  Array1::linspace(min_bl, max_bl, n_points).to_vec()
}
