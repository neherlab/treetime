use crate::make_report;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::params::BrentParams;
use crate::reroot::split::FindRootResult;
use crate::reroot::traits::RootStats;
use argmin::core::Executor;
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;

/// Optimize the split position along an edge with Brent's method on `[0, 1]`.
pub fn optimize_brent<S: RootStats>(
  edge: GraphEdgeKey,
  cost_fn: &EdgeCostFn<S>,
  params: &BrentParams,
) -> Result<FindRootResult<S>, Report> {
  let solver = BrentOpt::new(0.0, 1.0).set_tolerance(f64::EPSILON.sqrt(), params.brent_tolerance);

  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| cfg.max_iters(params.brent_max_iters as u64))
    .run()
    .map_err(|e| make_report!("Brent optimization failed: {e}"))?;

  let best_split = result
    .state
    .best_param
    .ok_or_else(|| make_report!("Brent optimization returned no parameter for edge {edge}"))?;
  let best_score = result.state.best_cost;
  let best_stats = cost_fn.evaluate(best_split);

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    stats: best_stats,
    score: best_score,
  })
}
