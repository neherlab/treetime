use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn leaf_order(graph: &Graph) -> Vec<GraphNodeKey> {
  graph.get_leaves().map(|leaf| leaf.key()).collect()
}
