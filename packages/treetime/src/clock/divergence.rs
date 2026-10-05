use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;

pub fn root_to_node_divergences(
  graph: &Graph,
  edge_length: impl Fn(GraphEdgeKey) -> f64,
) -> Result<BTreeMap<GraphNodeKey, f64>, Report> {
  let mut divergences = BTreeMap::new();
  graph.iter_depth_first_preorder_forward(|node| {
    let div = match node.parent_keys.first() {
      Some((parent_key, edge_key)) => divergences[parent_key] + edge_length(*edge_key),
      None => 0.0,
    };
    divergences.insert(node.key, div);
    Ok(())
  })?;
  Ok(divergences)
}

pub fn root_to_node_divergences_where_known(
  tree: &TreeView<'_>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> BTreeMap<GraphNodeKey, Option<f64>> {
  let mut divergences: BTreeMap<GraphNodeKey, Option<f64>> = BTreeMap::new();
  for &key in tree.preorder() {
    let divergence = match tree.parent(key) {
      Some((parent_key, edge_key)) => divergences[&parent_key]
        .zip(branch_lengths[&edge_key])
        .map(|(parent_divergence, length)| parent_divergence + length),
      None => Some(0.0),
    };
    divergences.insert(key, divergence);
  }
  divergences
}
