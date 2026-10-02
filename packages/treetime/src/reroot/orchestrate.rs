use crate::reroot::params::BrentParams;
use crate::reroot::placement::{RootTarget, leaf_keys, require_dated_new_leaves};
use crate::reroot::search::find_best_root;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{
  RerootResult, apply_reroot_topology, record_merge, record_split, remove_node_if_trivial, remove_stem_root,
  split_edge, trivial_node_branch_lengths,
};

pub(crate) fn reroot_in_place<S: RootStats>(
  graph: &mut Graph,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  root_stats: &S,
  variance: &VarianceModel,
  opt_params: &BrentParams,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<RerootResult, Report> {
  let best = find_best_root(graph, edge_stats, root_stats, variance, branch_lengths, opt_params)?;
  apply_root_at_edge(graph, best.edge, best.split, topo, branch_lengths, names)
}

pub(crate) fn reroot_at_node(
  graph: &mut Graph,
  node_key: GraphNodeKey,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<RerootResult, Report> {
  let edge = graph.parent_inbound_edge(node_key)?;
  apply_root_at_edge(graph, edge, 0.5, topo, branch_lengths, names)
}

fn apply_root_at_edge(
  graph: &mut Graph,
  edge: Option<GraphEdgeKey>,
  split: f64,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<RerootResult, Report> {
  let old_root_key = graph.root_key()?;
  let Some(edge_key) = edge else {
    return Ok(RerootResult::unchanged(old_root_key));
  };
  let target = RootTarget::on_edge(graph, edge_key, split, topo.split_edge)?;
  if target == RootTarget::Node(old_root_key) {
    return Ok(RerootResult::unchanged(old_root_key));
  }

  let leaves_before = leaf_keys(graph);
  let stem_removal = remove_stem_root(graph, old_root_key)?;
  let (current_root_key, target) = match &stem_removal {
    Some(stem) => {
      branch_lengths.remove(&stem.removed_edge_key);
      (stem.new_root_key, target.after_stem_removal(stem))
    },
    None => (old_root_key, target),
  };

  let (new_root_key, edge_split) = match target {
    RootTarget::Node(key) => (key, None),
    RootTarget::Split { edge_key, split } => {
      let length = branch_lengths.get(&edge_key).copied().flatten();
      let info = split_edge(graph, edge_key, split, length)?;
      record_split(branch_lengths, &info);
      (info.new_node_key, Some(info))
    },
  };

  let (inverted_edge_keys, edge_merge) = if new_root_key == current_root_key {
    (vec![], None)
  } else {
    let mut inverted = apply_reroot_topology(graph, current_root_key, new_root_key)?;
    let merge = if topo.remove_trivial_root {
      let (parent_branch, child_branch) = trivial_node_branch_lengths(graph, current_root_key, branch_lengths);
      remove_node_if_trivial(graph, current_root_key, parent_branch, child_branch)?
    } else {
      None
    };
    if let Some(merge) = &merge {
      inverted.retain(|k| *k != merge.parent_edge_key && *k != merge.child_edge_key);
      record_merge(branch_lengths, merge);
    }
    (inverted, merge)
  };

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
