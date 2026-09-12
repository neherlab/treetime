use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

/// The branch length that propagates sequence profiles along an edge: the clock-constrained length
/// when one has been committed, otherwise the raw ML or input length.
///
/// This is the domain choice that a timetree edge makes (`clock_branch_length` over the raw length);
/// every other command has no clock length and falls back to the raw length. Kept as a free function
/// so the choice stays named and testable where a profile map is derived.
pub fn profile_branch_length(clock: Option<f64>, raw: Option<f64>) -> Option<f64> {
  clock.or(raw)
}

/// Derive the per-edge profile branch length map (`f64`) each marginal pass propagates sequence
/// profiles along, from the raw input-tree branch length value map.
///
/// The value is `raw.unwrap_or(0.0)`: an edge with no length resolves to `0.0` for the passes. For a
/// timetree the clock-constrained length is combined in separately via
/// [`timetree_branch_lengths`](crate::timetree::inference::runner::timetree_branch_lengths); this
/// helper serves the non-timetree marginal passes (ancestral, optimize, clock, mugration).
pub fn profile_branch_lengths(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
  branch_lengths
    .iter()
    .map(|(key, raw)| (*key, profile_branch_length(None, *raw).unwrap_or(0.0)))
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
