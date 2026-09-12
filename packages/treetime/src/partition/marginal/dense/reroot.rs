use crate::ancestral::pipeline::DenseReconstruction;
use crate::partition::storage::dense::DenseNodeState;
use eyre::Report;
use treetime_graph::reroot::RerootChanges;

/// Apply a reroot to a completed dense reconstruction. The dense representation carries no per-edge
/// observation lists, so a reroot only changes the node set: it introduces the split node and removes
/// the merged trivial-root node. Seed a placeholder state for the introduced node, drop the removed
/// one, and clear the edge messages and estimates; the next marginal update rebuilds a complete set
/// over the rerooted topology from the leaf-seeded node states.
///
/// This is the value-style successor of the retired `PartitionRerootOps::apply_reroot` for dense.
pub fn reroot_dense(family: &mut DenseReconstruction, changes: &RerootChanges) -> Result<(), Report> {
  if let Some(info) = &changes.edge_split {
    family
      .node_states
      .entry(info.new_node_key)
      .or_insert_with(DenseNodeState::empty);
  }

  if let Some(info) = &changes.edge_merge {
    family.node_states.remove(&info.removed_node_key);
  }

  family.backward.clear();
  family.forward.clear();
  family.estimates.clear();

  Ok(())
}
