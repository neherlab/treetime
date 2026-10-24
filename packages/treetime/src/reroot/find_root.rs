use crate::clock::clock_regression::clock_regression;
use crate::clock::clock_state::ClockInputs;
use crate::progress::LogSink;
use crate::reroot::edge_search::search_edges;
use crate::reroot::midpoint::midpoint_root;
use crate::reroot::named_node::{named_node_point, oldest_tip, tip_group_mrca};
use crate::reroot::params::{RerootMethod, RerootParams, RerootSpec, RootObjective};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::split_fraction::SplitFraction;

pub(crate) fn find_root(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  reroot_params: &RerootParams,
  log: &dyn LogSink,
) -> Result<Option<EdgePoint>, Report> {
  let (force_positive, objective) = match &reroot_params.spec {
    RerootSpec::Method(RerootMethod::LeastSquares) => (reroot_params.force_positive_rate, reroot_params.objective),
    RerootSpec::Method(RerootMethod::MinDev) => (false, RootObjective::FixedRate(0.0)),
    RerootSpec::Method(RerootMethod::Oldest) => {
      return named_node_point(graph, branch_lengths, oldest_tip(graph, inputs)?);
    },
    RerootSpec::Tips(tips) => return named_node_point(graph, branch_lengths, tip_group_mrca(graph, tips)?),
    RerootSpec::Method(RerootMethod::Midpoint) => return midpoint_root(graph, branch_lengths),
  };
  let state = clock_regression(graph, inputs, outliers, branch_lengths, None)?;
  let best = search_edges(graph, &state, branch_lengths, force_positive, objective, log)?;
  Ok(best.point)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct EdgePoint {
  pub(crate) edge: GraphEdgeKey,
  pub(crate) split: SplitFraction,
}
