use crate::gtr::gtr::GTR;
use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::pipeline::SparseReconstruction;
use crate::make_internal_report;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{EdgeMergeInfo, RerootChanges};
use treetime_primitives::Seq;

/// Apply a reroot to the sparse inference state that survives it: rewrite the durable observations
/// (Fitch substitutions, indels, root sequence, and the split node's observations), then reconcile the
/// node-state map so the next marginal update sees exactly the current node set.
///
/// Takes the two values a reroot carries across and returns a reconstruction with no per-edge results:
/// the messages and estimates of the previous update describe the pre-reroot topology, and the next
/// marginal update rebuilds a complete set. The node states keep the leaf seeds and the derived root
/// sequence, so they are reconciled rather than dropped.
pub fn reroot_sparse(
  partition: PartitionMarginalSparse,
  gtr: GTR,
  node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  changes: &RerootChanges,
) -> Result<SparseReconstruction, Report> {
  let mut partition = partition;
  let mut node_states = node_states;

  apply_reroot_changes(&mut partition, &mut node_states, changes)?;

  if let Some(info) = &changes.edge_merge {
    remove_trivial_root(&mut partition, &mut node_states, info)?;
  }

  let node_states = reconcile_node_states(&partition, node_states);
  Ok(SparseReconstruction::seeded(partition, gtr, node_states))
}

// Phase 1: topology changes + root_sequence derivation + new node init.
// Must complete before remove_trivial_root (phase 2) because derive_root_sequence reads the
// parent_edge that phase 2 deletes.
//
// Note: root_composition + child-side fitch_subs + indels != child_composition. Non-char (N, gap)
// differences between nodes are not encoded as Fitch subs - they are tracked through non_char ranges
// on each node instead.
fn apply_reroot_changes(
  partition: &mut PartitionMarginalSparse,
  node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>,
  changes: &RerootChanges,
) -> Result<(), Report> {
  // Split edge: child-side gets all mutations, parent-side is empty
  if let Some(info) = &changes.edge_split {
    let old_edge_data = partition
      .obs_edges
      .remove(&info.old_edge_key)
      .ok_or_else(|| make_internal_report!("Old edge {:?} must exist for split", info.old_edge_key))?;
    partition.obs_edges.insert(info.child_side_edge_key, old_edge_data);
    partition
      .obs_edges
      .insert(info.parent_side_edge_key, SparseEdgeObs::default());
  }

  // Invert edges on the reroot path (substitutions and indels reverse direction).
  for edge_key in &changes.inverted_edge_keys {
    let edge_data = partition
      .obs_edges
      .get_mut(edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?;
    edge_data.invert_fitch_subs();
    for indel in &mut edge_data.indels {
      indel.invert();
    }
  }

  // Derive root_sequence for the new root
  derive_root_sequence(partition, changes);

  // Initialize the new split node's observations and node state from the finalized root_sequence.
  if let Some(info) = &changes.edge_split {
    partition.obs_nodes.insert(
      info.new_node_key,
      SparseNodeObs::new(&partition.root_sequence, &partition.alphabet),
    );
    node_states.insert(info.new_node_key, SparseNodeState::leaf(&partition.root_sequence));
  }

  Ok(())
}

fn remove_trivial_root(
  partition: &mut PartitionMarginalSparse,
  node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>,
  info: &EdgeMergeInfo,
) -> Result<(), Report> {
  let parent_edge = partition
    .obs_edges
    .remove(&info.parent_edge_key)
    .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
  let child_edge = partition
    .obs_edges
    .remove(&info.child_edge_key)
    .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

  partition.obs_nodes.remove(&info.removed_node_key);
  node_states.remove(&info.removed_node_key);

  let merged_subs = parent_edge.chain_fitch_subs(child_edge.fitch_subs())?;
  let merged_indels = parent_edge.chain_fitch_indels(&child_edge.indels);

  let mut merged_edge = SparseEdgeObs::default();
  merged_edge.set_fitch_subs(merged_subs);
  merged_edge.indels = merged_indels;

  partition.obs_edges.insert(info.merged_edge_key, merged_edge);
  Ok(())
}

fn derive_root_sequence(partition: &mut PartitionMarginalSparse, changes: &RerootChanges) {
  if !changes.inverted_edge_keys.is_empty() {
    let mut new_root_seq = partition.root_sequence.clone();
    for edge_key in &changes.inverted_edge_keys {
      if let Some(edge_data) = partition.obs_edges.get(edge_key) {
        apply_edge_to_sequence(&mut new_root_seq, edge_data, &partition.alphabet);
      }
    }
    partition.root_sequence = new_root_seq;
  } else if let Some(info) = &changes.edge_merge {
    if let Some(parent_edge) = partition.obs_edges.get(&info.parent_edge_key) {
      let parent_edge = parent_edge.clone();
      apply_edge_to_sequence(&mut partition.root_sequence, &parent_edge, &partition.alphabet);
    }
  }
}

fn apply_edge_to_sequence(seq: &mut Seq, edge: &SparseEdgeObs, alphabet: &Alphabet) {
  for sub in edge.fitch_subs() {
    if sub.pos() < seq.len() {
      seq[sub.pos()] = sub.reff();
    }
  }
  for indel in &edge.indels {
    if indel.range.0 < seq.len() && indel.range.1 <= seq.len() {
      if indel.is_deletion() {
        seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
      } else {
        seq[indel.range.0..indel.range.1].fill(alphabet.gap());
      }
    }
  }
}

/// Reconcile the node-state map to the partition's current node observations: seed a placeholder state
/// for every node introduced by the reroot and drop the state of any node it removed. Leaf seeds and the
/// derived root sequence are preserved; the marginal passes rebuild the evolving internal state.
fn reconcile_node_states(
  partition: &PartitionMarginalSparse,
  node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
) -> BTreeMap<GraphNodeKey, SparseNodeState> {
  let mut node_states = node_states;
  for key in partition.obs_nodes.keys().copied() {
    node_states.entry(key).or_insert_with(SparseNodeState::empty);
  }
  node_states.retain(|key, _| partition.obs_nodes.contains_key(key));
  node_states
}
