use crate::ancestral::pipeline::DenseReconstruction;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::storage::dense::DenseNodeState;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;

/// Apply a reroot to the dense inference state that survives it. The dense representation carries no
/// per-edge observation lists, so a reroot only changes the node set: it introduces the split node and
/// removes the merged trivial-root node.
///
/// Takes the two values a reroot carries across and returns a reconstruction with no per-edge results:
/// the messages and estimates of the previous update describe the pre-reroot topology, and the next
/// marginal update rebuilds a complete set from the leaf-seeded node states.
pub fn reroot_dense(
  partition: PartitionMarginalDense,
  node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  changes: &RerootChanges,
) -> DenseReconstruction {
  let mut node_states = node_states;

  if let Some(info) = &changes.edge_split {
    node_states
      .entry(info.new_node_key)
      .or_insert_with(DenseNodeState::empty);
  }

  if let Some(info) = &changes.edge_merge {
    node_states.remove(&info.removed_node_key);
  }

  DenseReconstruction::seeded(partition, node_states)
}
