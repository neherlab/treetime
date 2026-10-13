use crate::make_report;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::params::BrentParams;
use crate::reroot::split::FindRootResult;
use crate::reroot::traits::RootStats;
use argmin::core::Executor;
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_utils::cast::usize_to_u64;

pub(crate) fn optimize_brent<S: RootStats>(
  edge: GraphEdgeKey,
  cost_fn: &EdgeCostFn<S>,
  params: &BrentParams,
) -> Result<FindRootResult, Report> {
  let solver = BrentOpt::new(0.0, 1.0).set_tolerance(f64::EPSILON.sqrt(), params.brent_tolerance);

  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| cfg.max_iters(usize_to_u64(params.brent_max_iters)))
    .run()
    .map_err(|e| make_report!("Brent optimization failed: {e}"))?;

  let best_split = result
    .state
    .best_param
    .ok_or_else(|| make_report!("Brent optimization returned no parameter for edge {edge}"))?;
  let best_score = result.state.best_cost;

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    score: best_score,
  })
}
