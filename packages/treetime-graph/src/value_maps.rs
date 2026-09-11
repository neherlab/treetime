use crate::edge::GraphEdge;
use crate::graph::Graph;
use crate::node::{GraphNode, GraphNodeKey, Named};
use std::collections::BTreeMap;

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
