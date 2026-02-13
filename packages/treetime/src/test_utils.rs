use treetime_graph::edge::GraphEdge;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node::Named;

/// Find a node key by its name in a graph where nodes implement `Named`.
pub fn find_node_key_by_name<N, E, D>(graph: &Graph<N, E, D>, name: &str) -> Option<GraphNodeKey>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Send + Sync,
{
  for node in graph.get_nodes() {
    let node = node.read_arc();
    let payload = node.payload().read_arc();
    if payload.name().is_some_and(|n| n.as_ref() == name) {
      return Some(node.key());
    }
  }
  None
}

/// Find an edge key by source and target node names.
pub fn find_edge_key<N, E, D>(
  graph: &Graph<N, E, D>,
  source_name: &str,
  target_name: &str,
) -> Option<GraphEdgeKey>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: Send + Sync,
{
  let source_key = find_node_key_by_name(graph, source_name)?;
  let target_key = find_node_key_by_name(graph, target_name)?;

  for edge in graph.get_edges() {
    let edge = edge.read_arc();
    if edge.source() == source_key && edge.target() == target_key {
      return Some(edge.key());
    }
  }
  None
}
