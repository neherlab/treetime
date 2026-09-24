use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;

pub fn branch_lengths_or_zero(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
  branch_lengths
    .iter()
    .map(|(key, raw)| (*key, raw.unwrap_or(0.0)))
    .collect()
}

pub(crate) fn ancestral_reconstruction(
  graph: &Graph,
  mut advance: impl FnMut(&GraphNodeForward) -> Result<bool, Report>,
) -> Result<Vec<GraphNodeKey>, Report> {
  let mut emitted_nodes = Vec::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if advance(&node)? {
      emitted_nodes.push(node.key);
    }
    Ok(())
  })?;
  Ok(emitted_nodes)
}
