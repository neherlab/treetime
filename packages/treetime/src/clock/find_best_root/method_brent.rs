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

#[allow(
  clippy::as_conversions,
  clippy::unwrap_used,
  reason = "count/index numeric cast is exact for the domain range; unwrap on a value an upstream invariant guarantees is present"
)]
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
