use crate::clock::clock_regression::clock_regression;
use crate::clock::clock_state::{ClockEdgeInput, ClockInputs, ClockNodeInput, ClockState};
use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_root::FindRootResult;
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

/// What the root search fits and optimizes
pub(crate) struct RootSearch<'a> {
  pub outliers: &'a BTreeSet<GraphNodeKey>,
  pub reroot_params: &'a RerootParams,
  pub prev_clock_rate: Option<f64>,
}

/// Moves the root to the position that `search` selects. When the root moves and is an undated stem, the stem goes
/// first and the search repeats without it. The nodes and edges that the move creates get default clock inputs.
pub(crate) fn reroot_at_best_root(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &mut ClockInputs,
  search: &RootSearch<'_>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let mut best = search_root(graph, branch_lengths, inputs, search, log)?;
  let mut stem_removal = None;
  if root_moves(graph, best.edge, best.split, search.reroot_params.split_edge)? {
    stem_removal = remove_undated_stem(graph, branch_lengths, inputs)?;
    if stem_removal.is_some() {
      best = search_root(graph, branch_lengths, inputs, search, log)?;
    }
  }
  let result = apply_root_at_edge(
    graph,
    best.edge,
    best.split,
    search.reroot_params.topology(),
    branch_lengths,
    stem_removal,
    |key| inputs.likely_time(key).is_some(),
    names,
  )?;
  apply_reroot_to_inputs(inputs, &result);
  progress_info!(log, "Rerooted to {}", node_label(names, result.new_root_key));
  Ok(result)
}

/// Reroots a tree without dates by a `spec` that needs none (min-dev or tips). Every tip gets the same date and the
/// search runs at rate zero, so min-dev ranks the roots by the variance of the root-to-tip divergences.
pub(crate) fn reroot_undated(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  spec: &RerootSpec,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let tip_dates = graph.get_leaves().map(|leaf| (leaf.key(), Some(0.0))).collect();
  let mut inputs = ClockInputs::from_times(graph, &tip_dates, &BTreeMap::new());
  let reroot_params = RerootParams::new(spec.clone(), false).with_objective(RootObjective::FixedRate(0.0));
  let search = RootSearch {
    outliers: &BTreeSet::new(),
    reroot_params: &reroot_params,
    prev_clock_rate: None,
  };
  reroot_at_best_root(graph, branch_lengths, &mut inputs, &search, names, log)
}

fn search_root(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  inputs: &ClockInputs,
  search: &RootSearch<'_>,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  let state = clock_regression(graph, inputs, search.outliers, branch_lengths, search.prev_clock_rate)?;
  select_root(graph, inputs, &state, search.reroot_params, branch_lengths, log)
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

pub(crate) fn select_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  reroot_params: &RerootParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  match &reroot_params.spec {
    RerootSpec::Method(RerootMethod::LeastSquares | RerootMethod::ClockFilter) => find_best_root(
      graph,
      state,
      branch_lengths,
      reroot_params.force_positive_rate,
      reroot_params.objective,
      log,
    ),
    RerootSpec::Method(RerootMethod::MinDev) => {
      find_best_root(graph, state, branch_lengths, false, RootObjective::FixedRate(0.0), log)
    },
    RerootSpec::Method(RerootMethod::Oldest) => {
      find_oldest_root(graph, inputs, state, branch_lengths, reroot_params.objective)
    },
    RerootSpec::Tips(tips) => find_tip_group_root(graph, inputs, state, tips, branch_lengths, reroot_params.objective),
  }
}

#[derive(Clone, Debug, SmartDefault)]
pub struct RerootParams {
  pub(crate) spec: RerootSpec,

  pub(crate) objective: RootObjective,

  #[default = true]
  pub(crate) split_edge: bool,

  #[default = true]
  pub(crate) remove_trivial_root: bool,

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

  pub(crate) fn topology(&self) -> RerootTopologyParams {
    RerootTopologyParams {
      split_edge: self.split_edge,
      remove_trivial_root: self.remove_trivial_root,
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

fn find_oldest_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  let Some(oldest_key) = graph
    .get_leaves()
    .filter_map(|node| {
      let key = node.key();
      let time = inputs.likely_time(key)?;
      Some((time, key))
    })
    .min_by(|(lhs, _), (rhs, _)| lhs.total_cmp(rhs))
    .map(|(_, key)| key)
  else {
    return Err(input_error(
      "Cannot reroot to oldest tip because no dated leaves were found",
    ));
  };

  find_named_root_point(graph, inputs, state, oldest_key, branch_lengths, objective)
}

fn find_tip_group_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  tips: &[GraphNodeKey],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  if tips.is_empty() {
    return Err(input_error("--reroot-tips requires at least one tip name"));
  }

  let mrca_key = common_ancestor(graph, tips)?;
  find_named_root_point(graph, inputs, state, mrca_key, branch_lengths, objective)
}

fn find_named_root_point(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  node_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  let Some(edge) = graph.parent_inbound_edge(node_key)? else {
    let root = graph.get_exactly_one_root()?;
    let clock_set = state.node(root.key()).clone();
    return Ok(FindRootResult {
      edge: None,
      split: 0.0,
      chisq: objective.score(&clock_set),
      clock_set,
    });
  };

  let split = 0.5;
  let clock_set = BranchPointCostFunction::new(state, edge, branch_lengths, objective).evaluate_clock_set(split);
  Ok(FindRootResult {
    edge: Some(edge),
    split,
    chisq: objective.score(&clock_set),
    clock_set,
  })
}
