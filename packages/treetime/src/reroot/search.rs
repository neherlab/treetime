use crate::reroot::params::BrentParams;
use crate::reroot::split::{FindRootResult, find_best_split};
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Find the best root position over the whole tree.
///
/// Optimizes the split position on every edge and keeps the global minimum,
/// using the current root score as the baseline so an already-optimal tree is
/// left unchanged (`edge = None`). Optimizing each edge over `[0, 1]` covers
/// rooting at any existing node as a split endpoint, so no separate per-node
/// scan is needed.
pub fn find_best_root<N, E, D, S>(
  graph: &Graph<N, E, D>,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  root_stats: &S,
  variance: &VarianceModel,
  params: &BrentParams,
) -> Result<FindRootResult<S>, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
  S: RootStats,
{
  let mut best = FindRootResult {
    edge: None,
    split: 0.0,
    stats: root_stats.clone(),
    score: root_stats.score(),
  };

  for edge_obj in graph.get_edges() {
    let edge_key = edge_obj.read_arc().key();
    let res = find_best_split(graph, edge_key, edge_stats, variance, params)?;
    if res.score < best.score {
      best = res;
    }
  }

  Ok(best)
}
