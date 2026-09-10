use crate::edge::{GraphEdge, GraphEdgeKey, HasBranchLength};
use crate::graph::Graph;
use crate::node::{GraphNode, GraphNodeKey, Named};
use std::collections::BTreeMap;

/// Snapshot each edge's raw branch length into an edge-keyed value map.
///
/// The value is `HasBranchLength::branch_length()` verbatim, kept as `Option<f64>`: an edge with no
/// weight stays `None` rather than being coerced to a number, so a consumer that treats a missing
/// weight as an error keeps its error path. Distinct from
/// [`profile_branch_lengths`](../ancestral/marginal/fn.profile_branch_lengths.html), which resolves
/// to `f64` via the clock-constrained profile length for the marginal passes; this map is the raw
/// ML or input length that the ordering, GTR, divergence, clock, and writer readers consume.
pub fn edge_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> BTreeMap<GraphEdgeKey, Option<f64>>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .map(|edge| {
      let edge = edge.read_arc();
      (edge.key(), edge.payload().read_arc().branch_length())
    })
    .collect()
}

/// Snapshot each node's name into a node-keyed value map.
///
/// The value is `Named::name()` verbatim, kept as `Option<String>`: a node named only later by
/// [`assign_node_names`](../assign_node_names/fn.assign_node_names.html) stays `None` until it is
/// named, so a snapshot taken at a consumer's entry mirrors exactly what that consumer would read
/// off the payload at that point.
pub fn node_names<N, E, D>(graph: &Graph<N, E, D>) -> BTreeMap<GraphNodeKey, Option<String>>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Send + Sync,
{
  graph
    .get_nodes()
    .iter()
    .map(|node| {
      let node = node.read_arc();
      let name = node.payload().read_arc().name().map(|name| name.as_ref().to_owned());
      (node.key(), name)
    })
    .collect()
}
