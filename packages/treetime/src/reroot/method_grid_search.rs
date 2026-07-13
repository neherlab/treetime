use crate::make_error;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::params::GridSearchParams;
use crate::reroot::split::FindRootResult;
use crate::reroot::traits::RootStats;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;

/// Evaluate the objective on an equally-spaced grid of split positions and keep the best.
pub fn optimize_grid_search<S: RootStats>(
  edge: GraphEdgeKey,
  cost_fn: &EdgeCostFn<S>,
  params: &GridSearchParams,
) -> Result<FindRootResult<S>, Report> {
  let n_points = params.n_points;
  if n_points < 2 {
    return make_error!("Grid search requires at least 2 points, got {n_points}");
  }

  let mut best_score = f64::INFINITY;
  let mut best_split = f64::NAN;
  let mut best_stats = cost_fn.evaluate(0.0);

  // Equally-spaced split positions over `[0, 1]`, computed without allocating an
  // intermediate array. `step * i` reproduces `Array1::linspace(0.0, 1.0, n)`
  // bit-for-bit (linspace accumulates `start + step * i` with start = 0).
  let step = 1.0 / (n_points - 1) as f64;

  for i in 0..n_points {
    let x = step * i as f64;
    let stats = cost_fn.evaluate(x);
    let score = stats.score();
    if score < best_score {
      best_split = x;
      best_stats = stats;
      best_score = score;
    }
  }

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    stats: best_stats,
    score: best_score,
  })
}
