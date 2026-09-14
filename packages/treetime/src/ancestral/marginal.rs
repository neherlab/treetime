use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

/// Derive the per-edge branch length map (`f64`) a marginal pass propagates sequence profiles along,
/// from the raw input-tree branch length value map. An edge with no length resolves to `0.0`.
///
/// This is the raw-length collector: it applies no clock constraint. Once a timetree commit exists,
/// the timetree passes fold the committed clock length in instead via
/// [`timetree_branch_lengths`](crate::timetree::inference::runner::timetree_branch_lengths).
pub fn profile_branch_lengths(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
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
pub fn ancestral_reconstruction(
  graph: &Graph,
  mut reconstruct: impl FnMut(&GraphNodeForward) -> Option<Seq>,
  mut visitor: impl FnMut(GraphNodeKey, &Seq) -> Result<(), Report>,
) -> Result<BTreeMap<GraphNodeKey, Seq>, Report> {
  let mut node_sequences = BTreeMap::new();
  graph.iter_depth_first_preorder_forward(|node| match reconstruct(&node) {
    Some(seq) => {
      visitor(node.key, &seq)?;
      node_sequences.insert(node.key, seq);
      Ok(())
    },
    None => Ok(()),
  })?;
  Ok(node_sequences)
}
