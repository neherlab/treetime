use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::BrentParams;
use crate::make_report;
use crate::optimize::observer::OptimizationObserver;
use crate::progress::LogSink;
use crate::progress_info;
use argmin::core::Executor;
use argmin::core::observers::ObserverMode;
use argmin::solver::brent::BrentOpt;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_utils::cast::usize_to_u64;

pub(crate) fn optimize_brent(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &BrentParams,
  branch: &str,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  progress_info!(
    log,
    "Starting Brent optimization on the branch above {branch} with max_iters={}, tolerance={:.2e}",
    params.brent_max_iters,
    params.brent_tolerance
  );

  let solver = BrentOpt::new(0.0, 1.0);

  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| {
      cfg
        .max_iters(usize_to_u64(params.brent_max_iters))
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

  let best_split = result
    .state
    .best_param
    .ok_or_else(|| make_report!("Brent optimization returned no parameter"))?;
  let best_chisq = result.state.best_cost;

  progress_info!(
    log,
    "Brent optimization completed after {} iterations: best_split = {:.6}, best_cost = {:.6e}",
    result.state.iter,
    best_split,
    best_chisq
  );

  let best_clock_set = cost_fn.evaluate_clock_set(best_split)?;

  Ok(FindRootResult {
    edge: Some(edge),
    split: best_split,
    clock_set: best_clock_set,
    chisq: best_chisq,
  })
}
