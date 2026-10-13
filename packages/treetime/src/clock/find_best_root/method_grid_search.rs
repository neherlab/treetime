use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::GridSearchParams;
use crate::progress::LogSink;
use crate::progress_info;
use eyre::{Report, WrapErr};
use log::debug;
use ndarray::Array1;
use treetime_graph::edge::GraphEdgeKey;

pub(crate) fn optimize_grid_search(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &GridSearchParams,
  branch: &str,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  progress_info!(
    log,
    "Starting Grid Search optimization on the branch above {branch} with {} points",
    params.n_points
  );

  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_clock_set = cost_fn
    .evaluate_clock_set(0.0)
    .wrap_err("When evaluating the initial root position")?;

  let grid_points = Array1::linspace(0.0, 1.0, params.n_points);

  for (point_count, x) in grid_points.into_iter().enumerate() {
    let clock_set = cost_fn
      .evaluate_clock_set(x)
      .wrap_err_with(|| format!("When evaluating the root position {x}"))?;
    let chisq = cost_fn.score_clock_set(&clock_set);

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

  progress_info!(
    log,
    "Grid Search optimization completed: evaluated {} points, best_split = {:.6}, best_cost = {:.6e}",
    params.n_points,
    best_split,
    best_chisq
  );

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    clock_set: best_clock_set,
    chisq: best_chisq,
  })
}
