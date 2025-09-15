use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::BrentParams;
use crate::graph::edge::GraphEdgeKey;
use crate::make_report;
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{Executor, State};
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use log::{debug, info};

/// Observer for tracking Brent optimization progress
struct BrentObserver;

impl<I> Observe<I> for BrentObserver
where
  I: State,
  <I as State>::Param: std::fmt::Debug,
  <I as State>::Float: std::fmt::LowerExp,
{
  fn observe_iter(&mut self, state: &I, _kv: &argmin::core::KV) -> Result<(), argmin::core::Error> {
    // Only log every 10th iteration to reduce noise
    if state.get_iter() % 10 == 0 || state.get_iter() <= 5 {
      debug!(
        "Brent iteration {}: best_param = {:?}, best_cost = {:.6e}",
        state.get_iter(),
        state.get_best_param(),
        state.get_best_cost()
      );
    }
    Ok(())
  }
}

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
    .add_observer(BrentObserver, ObserverMode::Always)
    .run()
    .map_err(|e| make_report!("Brent optimization failed: {}", e))?;

  let best_split = result.state.best_param.unwrap();
  let best_chisq = result.state.best_cost;

  info!(
    "Brent optimization completed after {} iterations: best_split = {:.6}, best_cost = {:.6e}",
    result.state.iter, best_split, best_chisq
  );

  // Evaluate the cost function one more time to get the ClockSet data
  let best_total = cost_fn.evaluate_clock_set(best_split)?;

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    total: best_total,
    chisq: best_chisq,
  })
}
