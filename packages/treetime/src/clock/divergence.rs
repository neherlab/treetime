use eyre::Report;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::schedule::GraphSchedule;

pub fn root_to_node_divergences(
  graph: &Graph,
  edge_length: impl Fn(GraphEdgeKey) -> f64,
) -> Result<NodeMap<f64>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let mut divergences = NodeMap::new();
  for &key in schedule.preorder() {
    let div = match schedule.parents(key).first() {
      Some((parent_key, edge_key)) => divergences[parent_key] + edge_length(*edge_key),
      None => 0.0,
    };
    divergences.insert(key, div);
  }
  Ok(divergences)
}

pub fn root_to_node_divergences_where_known(
  schedule: &GraphSchedule<'_>,
  branch_lengths: &EdgeMap<Option<f64>>,
) -> Result<NodeMap<Option<f64>>, Report> {
  schedule.require_tree()?;
  let mut divergences: NodeMap<Option<f64>> = NodeMap::new();
  for &key in schedule.preorder() {
    let divergence = match schedule.parents(key).first() {
      Some((parent_key, edge_key)) => divergences[parent_key]
        .zip(branch_lengths[edge_key])
        .map(|(parent_divergence, length)| parent_divergence + length),
      None => Some(0.0),
    };
    divergences.insert(key, divergence);
  }
  Ok(divergences)
}
