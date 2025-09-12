use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::BrentParams;
use crate::graph::edge::GraphEdgeKey;
use crate::make_report;
use argmin::core::Executor;
use argmin::solver::brent::BrentOpt;
use eyre::Report;

/// Brent's method optimization for finding the best split point along an edge
pub fn optimize_brent(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &BrentParams,
) -> Result<FindRootResult, Report> {
  // Set up Brent solver with bounds [0.0, 1.0]
  // 0.0 means placing the root at the target node, 1.0 means placing it at the source node.
  let solver = BrentOpt::new(0.0, 1.0);

  // Run optimization
  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| {
      cfg
        .max_iters(params.brent_max_iters as u64)
        .target_cost(params.brent_tolerance)
    })
    .run()
    .map_err(|e| make_report!("Brent optimization failed: {}", e))?;

  let best_split = result.state.best_param.unwrap();
  let best_chisq = result.state.best_cost;

  // Evaluate the cost function one more time to get the ClockSet data
  let best_total = cost_fn.evaluate_clock_set(best_split)?;

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    total: best_total,
    chisq: best_chisq,
  })
}
