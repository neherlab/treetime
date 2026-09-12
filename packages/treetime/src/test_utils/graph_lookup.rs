use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub fn find_node_key_by_name<D>(
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  name: &str,
) -> Option<GraphNodeKey>
where
  D: Send + Sync,
{
  for node in graph.get_nodes() {
    let key = node.read_arc().key();
    if names.get(&key).and_then(|n| n.as_deref()) == Some(name) {
      return Some(key);
    }
  }
  None
}

pub fn find_edge_key<D>(
  graph: &Graph<D>,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  source_name: &str,
  target_name: &str,
) -> Option<GraphEdgeKey>
where
  D: Send + Sync,
{
  let source_key = find_node_key_by_name(graph, names, source_name)?;
  let target_key = find_node_key_by_name(graph, names, target_name)?;

  for edge in graph.get_edges() {
    let edge = edge.read_arc();
    if edge.source() == source_key && edge.target() == target_key {
      return Some(edge.key());
    }
  }
  None
}
