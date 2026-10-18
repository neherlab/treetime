use crate::clock::clock_regression::ClockTree;
use crate::clock::clock_state::{ClockEdgeInput, ClockInputs, ClockNodeInput, ClockState};
use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_root::FindRootResult;
use crate::clock::find_best_root::find_best_root::find_best_root;
use crate::clock::find_best_root::params::{RerootMethod, RerootSpec, RootObjective};
use crate::error::input_error;
use crate::progress::LogSink;
use crate::reroot::orchestrate::{RerootTopologyParams, apply_root_at_edge};
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{RerootResult, StemRemovalInfo, remove_stem_root};

/// Moves the root of `tree` to `best_root` and gives the nodes and edges that the move creates default clock inputs
pub(crate) fn reroot_clock_tree(
  mut tree: ClockTree,
  best_root: &FindRootResult,
  reroot_params: &RerootParams,
  stem_removal: Option<StemRemovalInfo>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(ClockTree, RerootResult), Report> {
  let ClockTree {
    graph,
    branch_lengths,
    inputs,
  } = &mut tree;
  let result = apply_root_at_edge(
    graph,
    best_root.edge,
    best_root.split,
    reroot_params.topology(),
    branch_lengths,
    stem_removal,
    |key| inputs.likely_time(key).is_some(),
    names,
  )?;
  apply_reroot_to_inputs(inputs, &result);
  Ok((tree, result))
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

pub(crate) fn remove_undated_stem(tree: &mut ClockTree) -> Result<Option<StemRemovalInfo>, Report> {
  let root_key = tree.graph.root_key()?;
  if tree.inputs.likely_time(root_key).is_some() {
    return Ok(None);
  }
  let Some(stem) = remove_stem_root(&mut tree.graph, root_key)? else {
    return Ok(None);
  };
  tree.branch_lengths.remove(&stem.removed_edge_key);
  tree.inputs.nodes.remove(&stem.removed_node_key);
  tree.inputs.edges.remove(&stem.removed_edge_key);
  Ok(Some(stem))
}

pub(crate) fn select_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  reroot_params: &RerootParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
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
