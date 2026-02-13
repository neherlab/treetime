use crate::commands::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::commands::clock::find_best_root::find_best_split::FindRootResult;
use crate::commands::clock::find_best_root::params::GoldenSectionParams;
use treetime_graph::edge::GraphEdgeKey;
use crate::make_report;
use argmin::core::observers::{Observe, ObserverMode};
use argmin::core::{Executor, State};
use argmin::solver::goldensectionsearch::GoldenSectionSearch;
use eyre::Report;
use log::{debug, info};

/// Observer for tracking Golden Section optimization progress
struct GoldenSectionObserver;

impl<I> Observe<I> for GoldenSectionObserver
where
  I: State,
  <I as State>::Param: std::fmt::Debug,
  <I as State>::Float: std::fmt::LowerExp,
{
  fn observe_iter(&mut self, state: &I, _kv: &argmin::core::KV) -> Result<(), argmin::core::Error> {
    // Only log every 10th iteration to reduce noise
    if state.get_iter().is_multiple_of(10) || state.get_iter() <= 5 {
      debug!(
        "Golden Section iteration {}: best_param = {:?}, best_cost = {:.6e}",
        state.get_iter(),
        state.get_best_param(),
        state.get_best_cost()
      );
    }
    Ok(())
  }
}

/// Golden section search optimization for finding the best split point along an edge
pub fn optimize_golden_section(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &GoldenSectionParams,
) -> Result<FindRootResult, Report> {
  info!(
    "Starting Golden Section optimization on edge {:?} with max_iters={}, tolerance={:.2e}",
    edge, params.golden_max_iters, params.golden_tolerance
  );
  // Set up Golden Section Search solver with bounds [0.0, 1.0]
  // 0.0 means placing the root at the target node, 1.0 means placing it at the source node.
  let solver = GoldenSectionSearch::new(0.0, 1.0)
    .map_err(|e| make_report!("Failed to create GoldenSectionSearch: {}", e))?
    .with_tolerance(params.golden_tolerance)
    .map_err(|e| make_report!("Golden Section optimization failed: {}", e))?;

  // Run optimization with initial guess at midpoint and observer
  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| {
      cfg
        .max_iters(params.golden_max_iters as u64)
        .target_cost(params.golden_tolerance)
        .param(0.5)
    })
    .add_observer(GoldenSectionObserver, ObserverMode::Always)
    .run()
    .map_err(|e| eyre::eyre!("Golden section search optimization failed: {}", e))?;

  let best_split = result.state.best_param.unwrap();
  let best_chisq = result.state.best_cost;

  info!(
    "Golden Section optimization completed after {} iterations: best_split = {:.6}, best_cost = {:.6e}",
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
