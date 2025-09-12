use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::GridSearchParams;
use crate::graph::edge::GraphEdgeKey;
use eyre::Report;
use ndarray::Array1;

/// Grid search optimization for finding the best split point along an edge
pub fn optimize_grid_search(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &GridSearchParams,
) -> Result<FindRootResult, Report> {
  // Grid search - interrogate different positions along the branch
  let mut best_chisq = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_total = cost_fn
    .evaluate_clock_set(0.0)
    .expect("Failed to evaluate initial position");

  for x in Array1::linspace(0.0, 1.0, params.n_points) {
    let clock_total = cost_fn.evaluate_clock_set(x).expect("Failed to evaluate clock set");
    let chisq = clock_total.chisq();
    if chisq < best_chisq {
      best_split = x;
      best_total = clock_total;
      best_chisq = chisq;
    }
  }

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    total: best_total,
    chisq: best_chisq,
  })
}
