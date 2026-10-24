use crate::clock::clock_regression::{ClockFit, ClockFitParams, ClockFitResult, ClockTree, fit_clock_model};
use crate::clock::clock_state::ClockInputs;
use crate::clock::date_constraints::DateConstraints;
use crate::progress::LogSink;
use crate::progress_info;
use crate::reroot::driver::reroot_at_best_root;
use crate::reroot::params::RerootParams;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::inference::result::given_times;
use eyre::{Report, WrapErr};
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::reroot::RerootResult;

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

  let DatedFit {
    tree: ClockTree {
      graph, branch_lengths, ..
    },
    fit,
    reroot,
  } = reroot_and_fit_to_dates(
    graph,
    branch_lengths,
    constraints,
    Some(params.reroot_params),
    &params.fit_params(),
    names,
    log,
  )
  .wrap_err("Failed to estimate clock model with reroot")?;

  let branch_model = match &reroot {
    Some(reroot) => branch_model.apply_reroot(&graph, &branch_lengths, reroot, log)?,
    None => branch_model,
  };

  Ok(RerootedTree {
    graph,
    branch_lengths,
    branch_model,
    clock_fit: fit.into_clock_fit()?,
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
      prev_clock_rate: self.prev_clock_rate,
    }
  }
}

pub(crate) fn reroot_and_fit_to_dates(
  graph: Graph,
  branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  constraints: &DateConstraints,
  reroot_params: Option<&RerootParams>,
  params: &ClockFitParams<'_>,
  names: &mut NodeNames,
  log: &dyn LogSink,
) -> Result<DatedFit, Report> {
  let times = given_times(&graph, constraints)?;
  let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
  let mut tree = ClockTree {
    graph,
    branch_lengths,
    inputs,
  };
  let reroot = match reroot_params {
    Some(reroot_params) => {
      let ClockTree {
        graph,
        branch_lengths,
        inputs,
      } = &mut tree;
      Some(reroot_at_best_root(
        graph,
        branch_lengths,
        inputs,
        params.outliers,
        &reroot_params.clone().with_clock_rate(params.clock_rate),
        names,
        log,
      )?)
    },
    None => None,
  };
  let fit = fit_clock_model(&tree, params, log)?;
  Ok(DatedFit { tree, fit, reroot })
}

pub(crate) struct DatedFit {
  pub tree: ClockTree,
  pub fit: ClockFitResult,
  pub reroot: Option<RerootResult>,
}

pub(crate) struct RerootedTree {
  pub graph: Graph,
  pub branch_lengths: BTreeMap<GraphEdgeKey, Option<f64>>,
  pub branch_model: BranchModel,
  pub clock_fit: ClockFit,
}
