use crate::reroot::div_stats::DivStats;
use crate::reroot::div_stats_traversal::compute_div_stats;
use crate::reroot::params::BrentParams;
use crate::reroot::placement::{RootTarget, leaf_keys, require_dated_new_leaves, root_moves};
use crate::reroot::search::find_best_root;
use crate::reroot::split::FindRootResult;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{
  RerootResult, StemRemovalInfo, apply_reroot_topology, record_merge, record_split, remove_node_if_trivial,
  remove_stem_root, split_edge, trivial_node_branch_lengths,
};

pub(crate) fn reroot_min_dev(
  graph: &mut Graph,
  variance: &VarianceModel,
  opt_params: &BrentParams,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<RerootResult, Report> {
  let best = min_dev_root(graph, variance, opt_params, branch_lengths)?;
  let stem_removal = if root_moves(graph, best.edge, best.split, topo.split_edge)? {
    remove_stem(graph, branch_lengths)?
  } else {
    None
  };
  let best = if stem_removal.is_some() {
    min_dev_root(graph, variance, opt_params, branch_lengths)?
  } else {
    best
  };
  apply_root_at_edge(graph, best.edge, best.split, topo, branch_lengths, stem_removal, names)
}

pub(crate) fn reroot_at_node(
  graph: &mut Graph,
  node_key: GraphNodeKey,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<RerootResult, Report> {
  let edge = graph.parent_inbound_edge(node_key)?;
  let stem_removal = if root_moves(graph, edge, 0.5, topo.split_edge)? {
    remove_stem(graph, branch_lengths)?
  } else {
    None
  };
  let edge = if stem_removal.is_some() {
    graph.parent_inbound_edge(node_key)?
  } else {
    edge
  };
  apply_root_at_edge(graph, edge, 0.5, topo, branch_lengths, stem_removal, names)
}

fn min_dev_root(
  graph: &Graph,
  variance: &VarianceModel,
  opt_params: &BrentParams,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<FindRootResult, Report> {
  let field = compute_div_stats(graph, branch_lengths, variance)?;
  find_best_root::<DivStats>(
    graph,
    &field.edge_stats,
    &field.root_stats,
    variance,
    branch_lengths,
    opt_params,
  )
}

fn remove_stem(
  graph: &mut Graph,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
) -> Result<Option<StemRemovalInfo>, Report> {
  let stem = remove_stem_root(graph, graph.root_key()?)?;
  if let Some(stem) = &stem {
    branch_lengths.remove(&stem.removed_edge_key);
  }
  Ok(stem)
}

fn apply_root_at_edge(
  graph: &mut Graph,
  edge: Option<GraphEdgeKey>,
  split: f64,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  stem_removal: Option<StemRemovalInfo>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<RerootResult, Report> {
  let old_root_key = graph.root_key()?;
  let unchanged = RerootResult {
    stem_removal: stem_removal.clone(),
    ..RerootResult::unchanged(old_root_key)
  };
  let Some(edge_key) = edge else {
    return Ok(unchanged);
  };
  let target = RootTarget::on_edge(graph, edge_key, split, topo.split_edge)?;
  if target == RootTarget::Node(old_root_key) {
    return Ok(unchanged);
  }

  let leaves_before = leaf_keys(graph);
  let (new_root_key, edge_split) = match target {
    RootTarget::Node(key) => (key, None),
    RootTarget::Split { edge_key, split } => {
      let length = branch_lengths.get(&edge_key).copied().flatten();
      let info = split_edge(graph, edge_key, split, length)?;
      record_split(branch_lengths, &info);
      (info.new_node_key, Some(info))
    },
  };

  let mut inverted_edge_keys = apply_reroot_topology(graph, old_root_key, new_root_key)?;
  let edge_merge = if topo.remove_trivial_root {
    let (parent_branch, child_branch) = trivial_node_branch_lengths(graph, old_root_key, branch_lengths);
    remove_node_if_trivial(graph, old_root_key, parent_branch, child_branch)?
  } else {
    None
  };
  if let Some(merge) = &edge_merge {
    inverted_edge_keys.retain(|k| *k != merge.parent_edge_key && *k != merge.child_edge_key);
    record_merge(branch_lengths, merge);
  }

  require_dated_new_leaves(graph, &leaves_before, |_| false, names)?;

  Ok(RerootResult {
    new_root_key,
    stem_removal,
    edge_split,
    edge_merge,
    inverted_edge_keys,
  })
}

#[derive(Debug, Clone, Copy, SmartDefault, Serialize, Deserialize)]
pub struct RerootTopologyParams {
  #[default = true]
  pub(crate) split_edge: bool,

  #[default = true]
  pub(crate) remove_trivial_root: bool,
}
