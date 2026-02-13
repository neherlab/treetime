use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::GridSearchParams;
use eyre::Report;
use log::{debug, info};
use ndarray::Array1;
use treetime_graph::edge::GraphEdgeKey;

/// Grid search optimization for finding the best split point along an edge
pub fn optimize_grid_search(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &GridSearchParams,
) -> Result<FindRootResult, Report> {
  info!(
    "Starting Grid Search optimization on edge {:?} with {} points",
    edge, params.n_points
  );

  // Grid search - interrogate different positions along the branch
  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_clock_set = cost_fn
    .evaluate_clock_set(0.0)
    .expect("Failed to evaluate initial position");

  let grid_points = Array1::linspace(0.0, 1.0, params.n_points);

  for (point_count, x) in grid_points.into_iter().enumerate() {
    let clock_set = cost_fn.evaluate_clock_set(x).expect("Failed to evaluate clock set");
    let chisq = clock_set.chisq();

    if chisq < best_chisq {
      debug!(
        "Grid search improved at point {}: x = {:.6}, chi-squared improved from {:.6e} to {:.6e}",
        point_count + 1,
        x,
        best_chisq,
        chisq
      );
      best_split = x;
      best_clock_set = clock_set;
      best_chisq = chisq;
    }
  }

  info!(
    "Grid Search optimization completed: evaluated {} points, best_split = {:.6}, best_cost = {:.6e}",
    params.n_points, best_split, best_chisq
  );

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    clock_set: best_clock_set,
    chisq: best_chisq,
  })
}
