use crate::clock::clock_regression::ClockVarianceParams;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{
  ClockEdgeInput, ClockEdgeState, ClockInputs, ClockNodeInput, ClockNodeState, ClockState,
};
use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_root::find_best_root;
use crate::clock::find_best_root::find_best_split::FindRootResult;
use crate::clock::find_best_root::params::{BranchPointOptimizationParams, RerootMethod, RerootSpec, RootObjective};
use crate::make_error;
use approx::ulps_eq;
use eyre::Report;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::common_ancestor::common_ancestor;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{
  EdgeSplitInfo, RerootResult, apply_reroot_topology, record_merge, record_split, remove_node_if_trivial, split_edge,
  trivial_node_branch_lengths,
};

#[derive(Clone, Debug, Serialize, Deserialize, SmartDefault)]
pub struct RerootParams {
  pub spec: RerootSpec,

  pub objective: RootObjective,

  #[default = true]
  pub split_edge: bool,

  #[default = true]
  pub remove_trivial_root: bool,

  #[default = true]
  pub force_positive_rate: bool,
}

impl RerootParams {
  #[must_use]
  pub(crate) fn with_objective(&self, objective: RootObjective) -> Self {
    Self {
      objective,
      ..self.clone()
    }
  }
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
pub(crate) fn reroot_in_place(
  graph: &mut Graph,
  inputs: &mut ClockInputs,
  mut state: ClockState,
  options: &ClockVarianceParams,
  params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(ClockState, RerootResult), Report> {
  let FindRootResult {
    edge, split, clock_set, ..
  } = select_root(
    graph,
    inputs,
    &state,
    options,
    params,
    reroot_params,
    branch_lengths,
    names,
  )?;

  let old_root_key = { graph.get_exactly_one_root()?.key() };
  let Some(edge_key) = edge else {
    return Ok((
      state,
      RerootResult {
        new_root_key: old_root_key,
        edge_split: None,
        edge_merge: None,
        inverted_edge_keys: vec![],
      },
    ));
  };

  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    (edge.source(), edge.target())
  };

  let (new_root_key, edge_split) = if ulps_eq!(split, 0.0, max_ulps = 5) {
    (source_key, None)
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    (target_key, None)
  } else if reroot_params.split_edge {
    let length = branch_lengths.get(&edge_key).copied().flatten();
    let split_info = create_new_root_node(graph, inputs, &mut state, edge_key, split, length, clock_set)?;
    record_split(branch_lengths, &split_info);
    (split_info.new_node_key, Some(split_info))
  } else {
    (if split < 0.5 { source_key } else { target_key }, None)
  };

  let (inverted_edge_keys, edge_merge) = if new_root_key != old_root_key {
    let mut inverted = apply_reroot(graph, &mut state, old_root_key, new_root_key, branch_lengths, options)?;

    let merge = if reroot_params.remove_trivial_root {
      let (parent_branch, child_branch) = trivial_node_branch_lengths(graph, old_root_key, branch_lengths);
      remove_node_if_trivial(graph, old_root_key, parent_branch, child_branch)?
    } else {
      None
    };

    if let Some(merge) = &merge {
      inverted.retain(|k| *k != merge.parent_edge_key && *k != merge.child_edge_key);
      record_merge(branch_lengths, merge);
      state.nodes.remove(&merge.removed_node_key);
      state.edges.remove(&merge.parent_edge_key);
      state.edges.remove(&merge.child_edge_key);
      state.edges.insert(merge.merged_edge_key, ClockEdgeState::default());
      inputs.nodes.remove(&merge.removed_node_key);
      inputs.edges.remove(&merge.parent_edge_key);
      inputs.edges.remove(&merge.child_edge_key);
      inputs.edges.insert(merge.merged_edge_key, ClockEdgeInput::default());
    }

    (inverted, merge)
  } else {
    (vec![], None)
  };

  Ok((
    state,
    RerootResult {
      new_root_key,
      edge_split,
      edge_merge,
      inverted_edge_keys,
    },
  ))
}

fn create_new_root_node(
  graph: &mut Graph,
  inputs: &mut ClockInputs,
  state: &mut ClockState,
  edge_key: GraphEdgeKey,
  split: f64,
  branch_length: Option<f64>,
  clock_set: ClockSet,
) -> Result<EdgeSplitInfo, Report> {
  let split_info = split_edge(graph, edge_key, split, branch_length)?;

  state.edges.remove(&split_info.old_edge_key);
  state
    .edges
    .insert(split_info.parent_side_edge_key, ClockEdgeState::default());
  state
    .edges
    .insert(split_info.child_side_edge_key, ClockEdgeState::default());
  state.nodes.insert(
    split_info.new_node_key,
    ClockNodeState {
      clock_set,
      ..ClockNodeState::default()
    },
  );

  inputs.edges.remove(&split_info.old_edge_key);
  inputs
    .edges
    .insert(split_info.parent_side_edge_key, ClockEdgeInput::default());
  inputs
    .edges
    .insert(split_info.child_side_edge_key, ClockEdgeInput::default());
  inputs.nodes.insert(split_info.new_node_key, ClockNodeInput::default());

  Ok(split_info)
}

fn select_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  options: &ClockVarianceParams,
  params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<FindRootResult, Report> {
  match &reroot_params.spec {
    RerootSpec::Method(RerootMethod::LeastSquares | RerootMethod::ClockFilter) => find_best_root(
      graph,
      inputs,
      state,
      options,
      params,
      branch_lengths,
      reroot_params.force_positive_rate,
      reroot_params.objective,
    ),
    RerootSpec::Method(RerootMethod::MinDev) => find_best_root(
      graph,
      inputs,
      state,
      options,
      params,
      branch_lengths,
      false,
      RootObjective::FixedRate(0.0),
    ),
    RerootSpec::Method(RerootMethod::Oldest) => {
      find_oldest_root(graph, inputs, state, options, branch_lengths, reroot_params.objective)
    },
    RerootSpec::Tips(tips) => find_tip_group_root(
      graph,
      inputs,
      state,
      options,
      tips,
      branch_lengths,
      reroot_params.objective,
      names,
    ),
  }
}

fn find_oldest_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  options: &ClockVarianceParams,
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
    return make_error!("Cannot reroot to oldest tip because no dated leaves were found");
  };

  find_named_root_point(graph, inputs, state, options, oldest_key, branch_lengths, objective)
}

fn find_tip_group_root(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  options: &ClockVarianceParams,
  tips: &[String],
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<FindRootResult, Report> {
  if tips.is_empty() {
    return make_error!("--reroot-tips requires at least one tip name");
  }

  let tip_keys = tips
    .iter()
    .map(|tip| {
      names
        .iter()
        .find(|(_, name)| name.as_deref() == Some(tip.as_str()))
        .map(|(key, _)| *key)
        .ok_or_else(|| eyre::eyre!("Reroot tip not found: {tip}"))
    })
    .try_collect::<_, Vec<_>, _>()?;
  let mrca_key = common_ancestor(graph, &tip_keys)?;
  find_named_root_point(graph, inputs, state, options, mrca_key, branch_lengths, objective)
}

fn find_named_root_point(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  options: &ClockVarianceParams,
  node_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  let Some(edge) = graph.parent_inbound_edge(node_key)? else {
    let root = graph.get_exactly_one_root()?;
    let clock_set = state.node(root.key()).clock_set.clone();
    return Ok(FindRootResult {
      edge: None,
      split: 0.0,
      chisq: objective.score(&clock_set),
      clock_set,
    });
  };

  let split = 0.5;
  let cost_fn = BranchPointCostFunction::new(graph, inputs, state, edge, branch_lengths, options, objective)?;
  let clock_set = cost_fn.evaluate_clock_set(split)?;
  Ok(FindRootResult {
    edge: Some(edge),
    split,
    chisq: objective.score(&clock_set),
    clock_set,
  })
}

#[allow(
  clippy::panic,
  clippy::unwrap_used,
  reason = "panics on a violated internal invariant; unwrap on a value an upstream invariant guarantees is present"
)]
fn apply_reroot(
  graph: &mut Graph,
  state: &mut ClockState,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &ClockVarianceParams,
) -> Result<Vec<GraphEdgeKey>, Report> {
  let inverted_edge_keys = apply_reroot_topology(graph, old_root_key, new_root_key)?;

  for edge_key in &inverted_edge_keys {
    let edge_len = branch_lengths[edge_key].unwrap();
    let branch_variance = options.variance_factor * edge_len + options.variance_offset;
    let edge_state = state
      .edges
      .get_mut(edge_key)
      .unwrap_or_else(|| panic!("Clock state is missing inverted edge {edge_key}"));
    let tmp_to_parent = edge_state.clock_to_parent.clone();
    edge_state.clock_to_parent = edge_state.clock_to_child.clone();
    edge_state.clock_to_child = tmp_to_parent;
    edge_state.clock_from_child = edge_state.clock_to_parent.propagate_averages(edge_len, branch_variance);
  }

  Ok(inverted_edge_keys)
}
