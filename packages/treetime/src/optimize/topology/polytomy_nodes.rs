use itertools::Itertools;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub fn find_polytomy_nodes(graph: &Graph) -> Vec<GraphNodeKey> {
  graph
    .get_nodes()
    .filter_map(|node| (node.degree_out() > 2).then_some(node.key()))
    .collect_vec()
}
