use crate::clock::clock_regression::{
  ClockFit, ClockFitParams, ClockRerootResult, ClockTree, estimate_clock_model_with_reroot_policy,
};
use crate::clock::clock_state::ClockInputs;
use crate::clock::date_constraints::DateConstraints;
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
use treetime_graph::node_names::NodeNames;

pub(crate) fn reroot_tree(
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  branch_model: BranchModel,
  constraints: &DateConstraints,
  params: &RerootFitParams<'_>,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<RerootedTree, Report> {
  progress_info!(
    log,
    "Reroot params: split_edge={}, remove_trivial_root={}, force_positive_rate={}",
    params.reroot_params.topology.split_edge,
    params.reroot_params.topology.remove_trivial_root,
    params.reroot_params.force_positive_rate
  );

  let (
    ClockTree {
      graph, branch_lengths, ..
    },
    clock_reroot_result,
  ) = fit_clock_to_dates(graph, branch_lengths, constraints, &params.fit_params(), names, log)
    .wrap_err("Failed to estimate clock model with reroot")?;

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

#[derive(Clone, Copy)]
pub(crate) struct RerootFitParams<'a> {
  pub outliers: &'a BTreeSet<GraphNodeKey>,
  pub clock_rate: Option<f64>,
  pub reroot_params: &'a RerootParams,
  pub prev_clock_rate: Option<f64>,
}

impl<'a> RerootFitParams<'a> {
  fn fit_params(&self) -> ClockFitParams<'a> {
    ClockFitParams {
      outliers: self.outliers,
      clock_rate: self.clock_rate,
      keep_root: false,
      reroot_params: self.reroot_params,
      prev_clock_rate: self.prev_clock_rate,
    }
  }
}

pub(crate) fn fit_clock_to_dates(
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  constraints: &DateConstraints,
  params: &ClockFitParams<'_>,
  names: &mut NodeNames,
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
    params,
    names,
    log,
  )
}

pub(crate) struct RerootedTree {
  pub graph: Graph,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub branch_model: BranchModel,
  pub clock_fit: ClockFit,
}
