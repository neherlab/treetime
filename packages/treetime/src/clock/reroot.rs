use crate::clock::clock_regression::ClockParams;
use crate::clock::clock_set::ClockSet;
use crate::clock::clock_state::{ClockEdgeInput, ClockEdgeState, ClockInputs, ClockNodeInput, ClockNodeState, ClockState};
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
  self as topology_reroot, record_merge, record_split, remove_node_if_trivial, split_edge, trivial_node_branch_lengths,
};

use topology_reroot::{EdgeSplitInfo, RerootResult};

/// Controls reroot behavior for graph topology changes.
#[derive(Clone, Debug, Serialize, Deserialize, SmartDefault)]
pub struct RerootParams {
  /// Root selection method.
  pub spec: RerootSpec,

  /// Objective used to rank candidate root positions for least-squares rerooting.
  pub objective: RootObjective,

  /// Allow creating a new node by splitting an edge during reroot.
  /// When false, reroot will snap to the nearest existing node endpoint.
  #[default = true]
  pub split_edge: bool,

  /// Remove the old root node if it becomes trivial (one parent, one child) after reroot.
  /// When false, the old root is preserved even if trivial.
  #[default = true]
  pub remove_trivial_root: bool,

  /// Only accept root positions with positive estimated clock rate.
  /// When false, the best chi-squared root is accepted regardless of rate sign.
  /// Use false for pre-filter steps where outliers may cause negative rates at all positions.
  #[default = true]
  pub force_positive_rate: bool,
}

impl RerootParams {
  #[must_use]
  pub fn with_objective(&self, objective: RootObjective) -> Self {
    Self {
      objective,
      ..self.clone()
    }
  }
}

pub fn reroot_in_place(
  graph: &mut Graph,
  inputs: &mut ClockInputs,
  mut state: ClockState,
  options: &ClockParams,
  params: &BranchPointOptimizationParams,
  reroot_params: &RerootParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(ClockState, RerootResult), Report> {
  let FindRootResult {
    edge, split, clock_set, ..
  } = select_root(graph, inputs, &state, options, params, reroot_params, branch_lengths, names)?;

  let old_root_key = { graph.get_exactly_one_root()?.read_arc().key() };
  let Some(edge_key) = edge else {
    // Already at the best root
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

  // Extract edge endpoints before the edge is removed by split_edge
  let (source_key, target_key) = {
    let edge = graph.get_edge(edge_key).expect("Edge not found");
    (edge.read_arc().source(), edge.read_arc().target())
  };

  // split = 0 roots at the source (parent), split = 1 at the target (child).
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

    // Remove edges consumed by the merge (they no longer exist in the graph)
    if let Some(merge) = &merge {
      inverted.retain(|k| *k != merge.parent_edge_key && *k != merge.child_edge_key);
      record_merge(branch_lengths, merge);
      // Keep the clock state and inputs consistent with the mutated node/edge set: the removed node
      // and the two edges it joined are gone; the merged edge takes their place with default messages
      // (recomputed by the next backward pass; a keep-root pass never reads them).
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

/// Create new root node by splitting the edge into two, then recording clock data for the new node
/// in the clock state.
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

  // Keep the clock state and inputs consistent with the mutated node/edge set: the split replaces one
  // edge with two and inserts one node. The new node carries the evaluated clock set at the split
  // point; its two edges start with default messages (recomputed by the next backward pass) and
  // default inputs.
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
  options: &ClockParams,
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
  options: &ClockParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  let Some(oldest_key) = graph
    .get_leaves()
    .into_iter()
    .filter_map(|node| {
      let key = node.read_arc().key();
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
  options: &ClockParams,
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
  options: &ClockParams,
  node_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  objective: RootObjective,
) -> Result<FindRootResult, Report> {
  let Some(edge) = graph.parent_inbound_edge(node_key)? else {
    let root = graph.get_exactly_one_root()?;
    let clock_set = state.node(root.read_arc().key()).clock_set.clone();
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

/// Modify graph topology to make the newly identified root the actual root,
/// then update clock-specific edge messages in the clock state.
fn apply_reroot(
  graph: &mut Graph,
  state: &mut ClockState,
  old_root_key: GraphNodeKey,
  new_root_key: GraphNodeKey,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  options: &ClockParams,
) -> Result<Vec<GraphEdgeKey>, Report> {
  let inverted_edge_keys = topology_reroot::apply_reroot_topology(graph, old_root_key, new_root_key)?;

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
