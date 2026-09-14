use itertools::Itertools;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

/// Find all nodes with more than 2 children (polytomies).
pub fn find_polytomy_nodes(graph: &Graph) -> Vec<GraphNodeKey> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node| (node.degree_out() > 2).then_some(node.key()))
    .collect_vec()
}
