use crate::make_report;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::method_brent::optimize_brent;
use crate::reroot::method_golden_section::optimize_golden_section;
use crate::reroot::method_grid_search::optimize_grid_search;
use crate::reroot::params::BranchPointOptimizationParams;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Outcome of a root search: the winning edge (or `None` for the current root),
/// the split fraction along it, the combined statistics, and the objective value.
#[derive(Debug, Clone)]
pub struct FindRootResult<S> {
  /// Edge carrying the best root position, or `None` when the current root wins.
  pub edge: Option<GraphEdgeKey>,

  /// Split fraction along the edge. `0` roots at the target node, `1` at the
  /// source node, `0 < x < 1` inserts a new node at that point.
  pub split: f64,

  /// Combined statistics at the winning position.
  pub stats: S,

  /// Objective value at the winning position (lower is better).
  pub score: f64,
}

/// Optimize the root position along a single edge using the configured 1D method.
pub fn find_best_split<N, E, D, S>(
  graph: &Graph<N, E, D>,
  edge: GraphEdgeKey,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  variance: &VarianceModel,
  params: &BranchPointOptimizationParams,
) -> Result<FindRootResult<S>, Report>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
  S: RootStats,
{
  let edge_obj = graph
    .get_edge(edge)
    .ok_or_else(|| make_report!("Edge not found: {edge}"))?;
  let (target_key, branch_length) = {
    let e = edge_obj.read_arc();
    let branch_length = e
      .payload()
      .read_arc()
      .branch_length()
      .ok_or_else(|| make_report!("Edge {edge} has no branch length"))?;
    (e.target(), branch_length)
  };

  let is_leaf = graph
    .get_node(target_key)
    .ok_or_else(|| make_report!("Target node not found for edge {edge}"))?
    .read_arc()
    .is_leaf();

  let (to_parent, to_child) = edge_stats
    .get(&edge)
    .ok_or_else(|| make_report!("No root statistics for edge {edge}"))?
    .clone();

  let cost_fn = EdgeCostFn {
    to_parent,
    to_child,
    branch_length,
    branch_variance: variance.branch(branch_length),
    is_leaf,
    leaf_time: None,
    variance_offset_leaf: variance.variance_offset_leaf,
  };

  match params {
    BranchPointOptimizationParams::Grid(p) => optimize_grid_search(edge, &cost_fn, p),
    BranchPointOptimizationParams::Brent(p) => optimize_brent(edge, &cost_fn, p),
    BranchPointOptimizationParams::GoldenSection(p) => optimize_golden_section(edge, &cost_fn, p),
  }
}
