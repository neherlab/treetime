use crate::seq::mutation::Sub;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;

pub fn compute_edge_mutation_counts(
  graph: &Graph,
  edge_subs: &BTreeMap<GraphEdgeKey, Vec<Sub>>,
) -> BTreeMap<GraphEdgeKey, usize> {
  graph
    .get_edges()
    .map(|edge| {
      let edge_key = edge.key();
      (edge_key, edge_subs[&edge_key].len())
    })
    .collect()
}
