use crate::clock::clock_regression::clock_regression;
use crate::clock::clock_state::ClockInputs;
use crate::clock::find_best_root::find_best_root::{FindRootResult, find_best_root};
use crate::clock::find_best_root::params::RootObjective;
use crate::progress::LogSink;
use crate::reroot::placement::{RootTarget, leaf_keys, require_dated_new_leaves, root_moves};
use eyre::Report;
use smart_default::SmartDefault;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{
  RerootResult, StemRemovalInfo, apply_reroot_topology, record_merge, record_split, remove_node_if_trivial,
  remove_stem_root, split_edge, trivial_node_branch_lengths,
};

/// Reroots at the position that minimizes the variance of the root-to-tip divergences
pub(crate) fn reroot_min_dev(
  graph: &mut Graph,
  topo: RerootTopologyParams,
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Result<RerootResult, Report> {
  let best = min_dev_root(graph, branch_lengths, log)?;
  let stem_removal = if root_moves(graph, best.edge, best.split, topo.split_edge)? {
    remove_stem(graph, branch_lengths)?
  } else {
    None
  };
  let best = if stem_removal.is_some() {
    min_dev_root(graph, branch_lengths, log)?
  } else {
    best
  };
  apply_root_at_edge(
    graph,
    best.edge,
    best.split,
    topo,
    branch_lengths,
    stem_removal,
    |_| false,
    names,
  )
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
  apply_root_at_edge(graph, edge, 0.5, topo, branch_lengths, stem_removal, |_| false, names)
}

/// The least-squares root at rate zero: with every tip at the same date, `chi^2 = (sum d^2 - (sum d)^2 / n) / 2n` is
/// the variance of the root-to-tip divergences over two
fn min_dev_root(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  log: &dyn LogSink,
) -> Result<FindRootResult, Report> {
  let tip_dates = graph.get_leaves().map(|leaf| (leaf.key(), Some(0.0))).collect();
  let inputs = ClockInputs::from_times(graph, &tip_dates, &BTreeMap::new());
  let state = clock_regression(graph, &inputs, &BTreeSet::new(), branch_lengths, None)?;
  find_best_root(graph, &state, branch_lengths, false, RootObjective::FixedRate(0.0), log)
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

/// Moves the root onto `edge` at fraction `split` from its parent, splitting the edge or snapping to an end as `topo`
/// says; a leaf that the move creates must satisfy `is_dated`
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
