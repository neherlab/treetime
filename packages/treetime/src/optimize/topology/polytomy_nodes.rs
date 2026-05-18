use itertools::Itertools;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey};

/// Find all nodes with more than 2 children (polytomies).
pub fn find_polytomy_nodes<N, E, D>(graph: &Graph<N, E, D>) -> Vec<GraphNodeKey>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node| {
      let node = node.read_arc();
      (node.degree_out() > 2).then_some(node.key())
    })
    .collect_vec()
}
