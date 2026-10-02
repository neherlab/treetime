use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::GoldenSectionParams;
use crate::make_report;
use crate::optimize::observer::OptimizationObserver;
use crate::progress::LogSink;
use crate::progress_info;
use argmin::core::Executor;
use argmin::core::observers::ObserverMode;
use argmin::solver::goldensectionsearch::GoldenSectionSearch;
use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;

#[allow(
  clippy::as_conversions,
  clippy::unwrap_used,
  reason = "count/index numeric cast is exact for the domain range; unwrap on a value an upstream invariant guarantees is present"
)]
pub(crate) fn optimize_golden_section(
  edge: GraphEdgeKey,
  cost_fn: &BranchPointCostFunction,
  params: &GoldenSectionParams,
  branch: &str,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  progress_info!(
    log,
    "Starting Golden Section optimization on the branch above {branch} with max_iters={}, tolerance={:.2e}",
    params.golden_max_iters,
    params.golden_tolerance
  );
  let solver = GoldenSectionSearch::new(0.0, 1.0)
    .map_err(|e| make_report!("Failed to create GoldenSectionSearch: {}", e))?
    .with_tolerance(params.golden_tolerance)
    .map_err(|e| make_report!("Golden Section optimization failed: {}", e))?;

  let result = Executor::new(cost_fn, solver)
    .configure(|cfg| {
      cfg
        .max_iters(params.golden_max_iters as u64)
        .target_cost(params.golden_tolerance)
        .param(0.5)
    })
    .add_observer(
      OptimizationObserver {
        label: "Golden Section",
        early_threshold: 5,
      },
      ObserverMode::Always,
    )
    .run()
    .map_err(|e| make_report!("Golden section search optimization failed: {e}"))?;

  let best_split = result.state.best_param.unwrap();
  let best_chisq = result.state.best_cost;

  progress_info!(
    log,
    "Golden Section optimization completed after {} iterations: best_split = {:.6}, best_cost = {:.6e}",
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
