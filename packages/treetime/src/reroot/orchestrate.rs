use crate::reroot::placement::{RootTarget, leaf_keys, require_dated_new_leaves};
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{
  RerootResult, StemRemovalInfo, apply_reroot_topology, record_merge, record_split, remove_node_if_trivial, split_edge,
  trivial_node_branch_lengths,
};

#[expect(
  clippy::too_many_arguments,
  reason = "each argument is an independent input of the move; a parameter struct would be built only for this call"
)]
pub(crate) fn apply_root_at_edge(
  graph: &mut Graph,
  edge: Option<GraphEdgeKey>,
  split: f64,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  stem_removal: Option<StemRemovalInfo>,
  is_dated: impl Fn(GraphNodeKey) -> bool,
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

  require_dated_new_leaves(graph, &leaves_before, is_dated, names)?;

  Ok(RerootResult {
    new_root_key,
    stem_removal,
    edge_split,
    edge_merge,
    inverted_edge_keys,
  })
}

#[derive(Debug, Clone, Copy, SmartDefault)]
pub struct RerootTopologyParams {
  #[default = true]
  pub(crate) split_edge: bool,

  #[default = true]
  pub(crate) remove_trivial_root: bool,
}
