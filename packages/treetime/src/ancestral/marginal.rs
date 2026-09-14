use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;

/// Resolve a per-edge branch length map to concrete `f64`, replacing a missing length with `0.0`.
pub fn branch_lengths_or_zero(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
  branch_lengths
    .iter()
    .map(|(key, raw)| (*key, raw.unwrap_or(0.0)))
    .collect()
}

/// Walk the graph in preorder, reconstructing each node's sequence with the supplied per-node
/// `reconstruct` closure, emitting each reconstructed sequence to `visitor`, and returning the
/// reconstructed sequences keyed by node id.
///
/// The representation-specific reconstruction (which reads the node states and, for sparse tips, the
/// forward messages) is supplied by the caller as a closure, so this single walk serves every marginal
/// representation without a shared mutable partition object. A node whose closure returns `None`
/// (a suppressed tip) is skipped.
/// Walk the tree in depth-first preorder, advancing each node's reconstruction state, and return the
/// node ids that emit a sequence, in walk order.
///
/// `advance` mutates the per-node state (recording posterior draws and tip reconstructions) and returns
/// `Some(())` when the node emits a sequence or `None` for a suppressed tip. The reconstructed sequences
/// themselves are read back off the partition afterward, so the walk holds no sequence in memory.
pub fn ancestral_reconstruction(
  graph: &Graph,
  mut advance: impl FnMut(&GraphNodeForward) -> Option<()>,
) -> Result<Vec<GraphNodeKey>, Report> {
  let mut emitted_nodes = Vec::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if advance(&node).is_some() {
      emitted_nodes.push(node.key);
    }
    Ok(())
  })?;
  Ok(emitted_nodes)
}
