use crate::make_report;
use crate::reroot::cost_function::EdgeCostFn;
use crate::reroot::method_brent::optimize_brent;
use crate::reroot::params::BrentParams;
use crate::reroot::traits::RootStats;
use crate::reroot::variance::VarianceModel;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

/// Outcome of a root search: the winning edge (or `None` for the current root),
/// the split fraction along it, the combined statistics, and the objective value.
#[derive(Debug, Clone)]
pub struct FindRootResult<S> {
  /// Edge carrying the best root position, or `None` when the current root wins.
  pub edge: Option<GraphEdgeKey>,

  /// Split fraction along the edge. `0` roots at the source (parent) node, `1`
  /// at the target (child) node, `0 < x < 1` inserts a new node at that point.
  pub split: f64,

  /// Combined statistics at the winning position.
  pub stats: S,

  /// Objective value at the winning position (lower is better).
  pub score: f64,
}

/// Optimize the root position along a single edge using Brent's method.
pub fn find_best_split<S>(
  graph: &Graph,
  edge: GraphEdgeKey,
  edge_stats: &BTreeMap<GraphEdgeKey, (S, S)>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  variance: &VarianceModel,
  params: &BrentParams,
) -> Result<FindRootResult<S>, Report>
where
  S: RootStats,
{
  let edge_obj = graph
    .get_edge(edge)
    .ok_or_else(|| make_report!("Edge not found: {edge}"))?;
  let target_key = edge_obj.read_arc().target();
  let branch_length = branch_lengths[&edge].ok_or_else(|| make_report!("Edge {edge} has no branch length"))?;

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

  optimize_brent(edge, &cost_fn, params)
}
