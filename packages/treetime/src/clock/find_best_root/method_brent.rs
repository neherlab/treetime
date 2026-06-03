use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::BrentParams;
use crate::make_report;
use crate::optimize::observer::OptimizationObserver;
use argmin::core::Executor;
use argmin::core::observers::ObserverMode;
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use log::info;
use treetime_graph::edge::GraphEdgeKey;

/// Brent's method optimization for finding the best split point along an edge
pub fn optimize_brent(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &BrentParams,
) -> Result<FindRootResult, Report> {
  info!(
    "Starting Brent optimization on edge {:?} with max_iters={}, tolerance={:.2e}",
    edge, params.brent_max_iters, params.brent_tolerance
  );

  // Set up Brent solver with bounds [0.0, 1.0]
  // 0.0 means placing the root at the target node, 1.0 means placing it at the source node.
  let solver = BrentOpt::new(0.0, 1.0);

  // Run optimization with observer
  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| {
      cfg
        .max_iters(params.brent_max_iters as u64)
        .target_cost(params.brent_tolerance)
    })
    .add_observer(
      OptimizationObserver {
        label: "Brent",
        early_threshold: 5,
      },
      ObserverMode::Always,
    )
    .run()
    .map_err(|e| make_report!("Brent optimization failed: {}", e))?;

  let best_split = result.state.best_param.unwrap();
  let best_chisq = result.state.best_cost;

  info!(
    "Brent optimization completed after {} iterations: best_split = {:.6}, best_cost = {:.6e}",
    result.state.iter, best_split, best_chisq
  );

  // Evaluate the cost function one more time to get the ClockSet data
  let best_clock_set = cost_fn.evaluate_clock_set(best_split)?;

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    clock_set: best_clock_set,
    chisq: best_chisq,
  })
}
