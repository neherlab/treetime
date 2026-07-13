use crate::make_report;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::params::GoldenSectionParams;
use crate::reroot::split::FindRootResult;
use crate::reroot::traits::RootStats;
use argmin::core::Executor;
use argmin::solver::goldensectionsearch::GoldenSectionSearch;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;

/// Optimize the split position along an edge with golden-section search on `[0, 1]`.
pub fn optimize_golden_section<S: RootStats>(
  edge: GraphEdgeKey,
  cost_fn: &EdgeCostFn<S>,
  params: &GoldenSectionParams,
) -> Result<FindRootResult<S>, Report> {
  let solver = GoldenSectionSearch::new(0.0, 1.0)
    .map_err(|e| make_report!("Failed to create GoldenSectionSearch: {e}"))?
    .with_tolerance(params.golden_tolerance)
    .map_err(|e| make_report!("Failed to configure GoldenSectionSearch: {e}"))?;

  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| {
      cfg
        .max_iters(params.golden_max_iters as u64)
        .target_cost(params.golden_tolerance)
        .param(0.5)
    })
    .run()
    .map_err(|e| make_report!("Golden section optimization failed: {e}"))?;

  let best_split = result
    .state
    .best_param
    .ok_or_else(|| make_report!("Golden section optimization returned no parameter for edge {edge}"))?;
  let best_score = result.state.best_cost;
  let best_stats = cost_fn.evaluate(best_split);

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    stats: best_stats,
    score: best_score,
  })
}
