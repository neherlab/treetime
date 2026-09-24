use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn find_node_key_by_name(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  name: &str,
) -> Option<GraphNodeKey> {
  for node in graph.get_nodes() {
    let key = node.key();
    if names.get(&key).and_then(|n| n.as_deref()) == Some(name) {
      return Some(key);
    }
  }
  None
}

pub(crate) fn find_edge_key(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  source_name: &str,
  target_name: &str,
) -> Option<GraphEdgeKey> {
  let source_key = find_node_key_by_name(graph, names, source_name)?;
  let target_key = find_node_key_by_name(graph, names, target_name)?;

  for edge in graph.get_edges() {
    if edge.source() == source_key && edge.target() == target_key {
      return Some(edge.key());
    }
  }
  None
}
