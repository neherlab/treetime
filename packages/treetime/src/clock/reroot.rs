use crate::clock::clock_regression::clock_regression;
use crate::clock::clock_state::{ClockEdgeInput, ClockInputs, ClockNodeInput};
use crate::clock::find_best_root::find_best_root::find_best_root;
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec, RootObjective};
use crate::error::input_error;
use crate::node_label::node_label;
use crate::progress::LogSink;
use crate::progress_info;
use crate::reroot::orchestrate::{RerootTopologyParams, apply_root_at_edge};
use crate::reroot::placement::root_moves;
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{RerootResult, StemRemovalInfo, remove_stem_root};

pub(crate) fn reroot_undated(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  spec: &RerootSpec,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let tip_dates = graph.get_leaves().map(|leaf| (leaf.key(), Some(0.0))).collect();
  let mut inputs = ClockInputs::from_times(graph, &tip_dates, &BTreeMap::new());
  let reroot_params = RerootParams::new(spec.clone(), false);
  reroot_at_best_root(
    graph,
    branch_lengths,
    &mut inputs,
    &BTreeSet::new(),
    &reroot_params,
    names,
    log,
  )
}

pub(crate) fn reroot_at_best_root(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &mut ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  reroot_params: &RerootParams,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let split_edge = reroot_params.topology.split_edge;
  let (mut edge, mut split) = select_root(graph, branch_lengths, inputs, outliers, reroot_params, log)?;
  let mut stem_removal = None;
  if root_moves(graph, edge, split, split_edge)? {
    stem_removal = remove_undated_stem(graph, branch_lengths, inputs)?;
    if stem_removal.is_some() {
      (edge, split) = select_root(graph, branch_lengths, inputs, outliers, reroot_params, log)?;
    }
  }
  let result = apply_root_at_edge(
    graph,
    edge,
    split,
    reroot_params.topology,
    branch_lengths,
    stem_removal,
    |key| inputs.likely_time(key).is_some(),
    names,
  )?;
  apply_reroot_to_inputs(inputs, &result);
  progress_info!(log, "Rerooted to {}", node_label(names, result.new_root_key));
  Ok(result)
}

fn apply_reroot_to_inputs(inputs: &mut ClockInputs, result: &RerootResult) {
  if let Some(split) = &result.edge_split {
    inputs.edges.remove(&split.old_edge_key);
    inputs
      .edges
      .insert(split.parent_side_edge_key, ClockEdgeInput::default());
    inputs
      .edges
      .insert(split.child_side_edge_key, ClockEdgeInput::default());
    inputs.nodes.insert(split.new_node_key, ClockNodeInput::default());
  }
  if let Some(merge) = &result.edge_merge {
    inputs.nodes.remove(&merge.removed_node_key);
    inputs.edges.remove(&merge.parent_edge_key);
    inputs.edges.remove(&merge.child_edge_key);
    inputs.edges.insert(merge.merged_edge_key, ClockEdgeInput::default());
  }
}

fn remove_undated_stem(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &mut ClockInputs,
) -> Result<Option<StemRemovalInfo>, Report> {
  let root_key = graph.root_key()?;
  if inputs.likely_time(root_key).is_some() {
    return Ok(None);
  }
  let Some(stem) = remove_stem_root(graph, root_key)? else {
    return Ok(None);
  };
  branch_lengths.remove(&stem.removed_edge_key);
  inputs.nodes.remove(&stem.removed_node_key);
  inputs.edges.remove(&stem.removed_edge_key);
  Ok(Some(stem))
}

fn select_root(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &ClockInputs,
  outliers: &BTreeSet<GraphNodeKey>,
  reroot_params: &RerootParams,
  log: &dyn LogSink,
) -> Result<(Option<GraphEdgeKey>, f64), Report> {
  let (force_positive, objective) = match &reroot_params.spec {
    RerootSpec::Method(RerootMethod::LeastSquares | RerootMethod::ClockFilter) => {
      (reroot_params.force_positive_rate, reroot_params.objective)
    },
    RerootSpec::Method(RerootMethod::MinDev) => (false, RootObjective::FixedRate(0.0)),
    RerootSpec::Method(RerootMethod::Oldest) => return named_root_point(graph, oldest_tip(graph, inputs)?),
    RerootSpec::Tips(tips) => return named_root_point(graph, tip_group_mrca(graph, tips)?),
  };
  let state = clock_regression(graph, inputs, outliers, branch_lengths, None)?;
  let best = find_best_root(graph, &state, branch_lengths, force_positive, objective, log)?;
  Ok((best.edge, best.split))
}

#[derive(Clone, Debug, SmartDefault)]
pub struct RerootParams {
  pub(crate) spec: RerootSpec,

  pub(crate) objective: RootObjective,

  pub(crate) topology: RerootTopologyParams,

  #[default = true]
  pub(crate) force_positive_rate: bool,
}

impl RerootParams {
  #[must_use]
  pub(crate) fn new(spec: RerootSpec, force_positive_rate: bool) -> Self {
    Self {
      spec,
      force_positive_rate,
      ..Self::default()
    }
  }

  #[must_use]
  pub(crate) fn with_objective(&self, objective: RootObjective) -> Self {
    Self {
      objective,
      ..self.clone()
    }
  }
}

fn oldest_tip(graph: &Graph, inputs: &ClockInputs) -> Result<GraphNodeKey, Report> {
  graph
    .get_leaves()
    .filter_map(|node| {
      let key = node.key();
      let time = inputs.likely_time(key)?;
      Some((time, key))
    })
    .min_by(|(lhs, _), (rhs, _)| lhs.total_cmp(rhs))
    .map(|(_, key)| key)
    .ok_or_else(|| input_error("Cannot reroot to oldest tip because no dated leaves were found"))
}

fn tip_group_mrca(graph: &Graph, tips: &[GraphNodeKey]) -> Result<GraphNodeKey, Report> {
  if tips.is_empty() {
    return Err(input_error("--reroot-tips requires at least one tip name"));
  }
  common_ancestor(graph, tips)
}

fn named_root_point(graph: &Graph, node_key: GraphNodeKey) -> Result<(Option<GraphEdgeKey>, f64), Report> {
  Ok(match graph.parent_inbound_edge(node_key)? {
    Some(edge) => (Some(edge), 0.5),
    None => (None, 0.0),
  })
}
