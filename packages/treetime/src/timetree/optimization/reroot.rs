use crate::clock::clock_regression::{
  ClockFit, ClockRerootResult, ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy,
};
use crate::clock::clock_state::ClockInputs;
use crate::clock::date_constraints::DateConstraints;
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::progress::LogSink;
use crate::progress_info;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::inference::result::given_times;
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of this step; a parameter struct would be built only for this call"
)]
pub(crate) fn reroot_tree(
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  branch_model: BranchModel,
  constraints: &DateConstraints,
  outliers: &BTreeSet<GraphNodeKey>,
  clock_params: &ClockVarianceParams,
  clock_rate: Option<f64>,
  branch_params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<RerootedTree, Report> {
  progress_info!(
    log,
    "Reroot params: split_edge={}, remove_trivial_root={}, force_positive_rate={}",
    reroot_params.split_edge,
    reroot_params.remove_trivial_root,
    reroot_params.force_positive_rate
  );

  let fit = DatedClockFit {
    outliers,
    clock_params,
    clock_rate,
    keep_root: false,
    branch_params,
    reroot_params,
    failure: "Failed to estimate clock model with reroot",
  };
  let (
    ClockTree {
      graph, branch_lengths, ..
    },
    clock_reroot_result,
  ) = fit_clock_to_dates(graph, branch_lengths, constraints, &fit, names, log)?;

  let branch_model = match clock_reroot_result.reroot_result() {
    Some(reroot) => branch_model.apply_reroot(&graph, &branch_lengths, reroot, log)?,
    None => branch_model,
  };

  Ok(RerootedTree {
    graph,
    branch_lengths,
    branch_model,
    clock_fit: clock_reroot_result.into_clock_fit()?,
  })
}

pub(crate) fn fit_clock_to_dates(
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  constraints: &DateConstraints,
  fit: &DatedClockFit<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<(ClockTree, ClockRerootResult), Report> {
  let times = given_times(&graph, constraints)?;
  let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
  estimate_clock_model_with_reroot_policy(
    ClockTree {
      graph,
      branch_lengths,
      inputs,
    },
    fit.outliers,
    fit.clock_params,
    fit.clock_rate,
    fit.keep_root,
    fit.branch_params,
    fit.reroot_params,
    None,
    names,
    log,
  )
  .wrap_err(fit.failure)
}

pub(crate) struct DatedClockFit<'a> {
  pub outliers: &'a BTreeSet<GraphNodeKey>,
  pub clock_params: &'a ClockVarianceParams,
  pub clock_rate: Option<f64>,
  pub keep_root: bool,
  pub branch_params: &'a BranchPointOptimizationParams,
  pub reroot_params: &'a RerootParams,
  pub failure: &'static str,
}

pub(crate) struct RerootedTree {
  pub graph: Graph,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub branch_model: BranchModel,
  pub clock_fit: ClockFit,
}
