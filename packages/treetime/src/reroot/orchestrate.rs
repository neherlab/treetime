use crate::make_report;
use crate::reroot::params::BrentParams;
use crate::reroot::search::find_best_root;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use approx::ulps_eq;
use eyre::Report;
use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};
use treetime_graph::reroot::{RerootResult, apply_reroot_topology, remove_node_if_trivial, split_edge};

/// Topology options applied once a root position has been chosen.
#[derive(Debug, Clone, Copy, SmartDefault, Serialize, Deserialize)]
pub struct RerootTopologyParams {
  /// Insert a new node when the chosen root falls in the interior of an edge.
  /// When false, snap to the nearer endpoint instead.
  #[default = true]
  pub split_edge: bool,

  /// Remove the old root if it becomes a trivial degree-2 node after rerooting.
  #[default = true]
  pub remove_trivial_root: bool,
}

/// Search for the best root by `S` scoring and apply it to the graph topology.
///
/// `fixup` runs after edges are inverted, receiving the inverted edge keys, to
/// repair domain-specific edge data. Objectives whose statistics are ephemeral
/// (e.g. divergence-only rooting) pass a no-op; partition state is reconciled
/// separately by the caller from the returned `RerootResult`.
pub fn reroot_in_place<N, E, D, S, F>(
  graph: &mut Graph<N, E, D>,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  root_stats: &S,
  variance: &VarianceModel,
  opt_params: &BrentParams,
  topo: RerootTopologyParams,
  fixup: F,
) -> Result<RerootResult, Report>
where
  N: GraphNode + Default,
  E: GraphEdge + HasBranchLength + Default,
  D: Send + Sync,
  S: RootStats,
  F: FnMut(&mut Graph<N, E, D>, &[GraphEdgeKey]) -> Result<(), Report>,
{
  let best = find_best_root(graph, edge_stats, root_stats, variance, opt_params)?;
  apply_root_at_edge(graph, best.edge, best.split, topo, fixup)
}

/// Reroot on the branch leading to `node_key` at its midpoint.
///
/// Used for tip- or MRCA-based rerooting, which needs no scoring. When the node
/// is already the root, the tree is left unchanged.
pub fn reroot_at_node<N, E, D, F>(
  graph: &mut Graph<N, E, D>,
  node_key: GraphNodeKey,
  topo: RerootTopologyParams,
  fixup: F,
) -> Result<RerootResult, Report>
where
  N: GraphNode + Default,
  E: GraphEdge + HasBranchLength + Default,
  D: Send + Sync,
  F: FnMut(&mut Graph<N, E, D>, &[GraphEdgeKey]) -> Result<(), Report>,
{
  let edge = graph.parent_inbound_edge(node_key)?;
  apply_root_at_edge(graph, edge, 0.5, topo, fixup)
}

fn apply_root_at_edge<N, E, D, F>(
  graph: &mut Graph<N, E, D>,
  edge: Option<GraphEdgeKey>,
  split: f64,
  topo: RerootTopologyParams,
  mut fixup: F,
) -> Result<RerootResult, Report>
where
  N: GraphNode + Default,
  E: GraphEdge + HasBranchLength + Default,
  D: Send + Sync,
  F: FnMut(&mut Graph<N, E, D>, &[GraphEdgeKey]) -> Result<(), Report>,
{
  let old_root_key = graph.get_exactly_one_root()?.read_arc().key();

  let Some(edge_key) = edge else {
    return Ok(RerootResult {
      new_root_key: old_root_key,
      edge_split: None,
      edge_merge: None,
      inverted_edge_keys: vec![],
    });
  };

  let (source_key, target_key) = {
    let edge = graph
      .get_edge(edge_key)
      .ok_or_else(|| make_report!("Edge not found: {edge_key}"))?;
    let edge = edge.read_arc();
    (edge.source(), edge.target())
  };

  // split = 0 roots at the source (parent), split = 1 at the target (child).
  let (new_root_key, edge_split) = if ulps_eq!(split, 0.0, max_ulps = 5) {
    (source_key, None)
  } else if ulps_eq!(split, 1.0, max_ulps = 5) {
    (target_key, None)
  } else if topo.split_edge {
    let info = split_edge(graph, edge_key, split)?;
    (info.new_node_key, Some(info))
  } else {
    (if split < 0.5 { source_key } else { target_key }, None)
  };

  let (inverted_edge_keys, edge_merge) = if new_root_key != old_root_key {
    let mut inverted = apply_reroot_topology(graph, old_root_key, new_root_key)?;
    fixup(graph, &inverted)?;

    let merge = if topo.remove_trivial_root {
      remove_node_if_trivial(graph, old_root_key)?
    } else {
      None
    };

    if let Some(merge) = &merge {
      inverted.retain(|k| *k != merge.parent_edge_key && *k != merge.child_edge_key);
    }

    (inverted, merge)
  } else {
    (vec![], None)
  };

  Ok(RerootResult {
    new_root_key,
    edge_split,
    edge_merge,
    inverted_edge_keys,
  })
}
