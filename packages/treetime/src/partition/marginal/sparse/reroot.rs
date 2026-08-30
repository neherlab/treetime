use crate::alphabet::alphabet::Alphabet;
use crate::make_internal_report;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition};
use crate::partition::traits::PartitionRerootOps;
use eyre::Report;
use treetime_graph::reroot::{EdgeMergeInfo, RerootChanges};
use treetime_primitives::Seq;

// Split by codepath stage: these helpers sit beside the trait impl that consumes them.
#[allow(clippy::multiple_inherent_impl)]
impl PartitionMarginalSparse {
  // Phase 1: topology changes + root_sequence derivation + new node init.
  // Must complete before remove_trivial_root (phase 2) because
  // derive_root_sequence reads the parent_edge that phase 2 deletes.
  //
  // Note: root_composition + child-side fitch_subs + indels != child_composition.
  // Non-char (N, gap) differences between nodes are not encoded as Fitch subs -
  // they are tracked through non_char ranges on each node instead.
  fn apply_reroot_changes(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    // Split edge: child-side gets all mutations, parent-side is empty
    if let Some(info) = &changes.edge_split {
      let mut old_edge_data = self
        .edges
        .remove(&info.old_edge_key)
        .ok_or_else(|| make_internal_report!("Old edge {:?} must exist for split", info.old_edge_key))?;
      old_edge_data.clear_ml_subs();
      self.edges.insert(info.child_side_edge_key, old_edge_data);
      self
        .edges
        .insert(info.parent_side_edge_key, SparseEdgePartition::default());
    }

    // Invert edges on the reroot path
    for edge_key in &changes.inverted_edge_keys {
      let edge_data = self
        .edges
        .get_mut(edge_key)
        .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?;
      edge_data.invert_for_reroot();
    }

    // Derive root_sequence for the new root
    self.derive_root_sequence(changes);

    // Initialize the new split node from the finalized root_sequence
    if let Some(info) = &changes.edge_split {
      self.nodes.insert(
        info.new_node_key,
        SparseNodePartition::new(&self.root_sequence, &self.alphabet)?,
      );
    }

    Ok(())
  }

  fn remove_trivial_root(&mut self, info: &EdgeMergeInfo) -> Result<(), Report> {
    let parent_edge = self
      .edges
      .remove(&info.parent_edge_key)
      .ok_or_else(|| make_internal_report!("Parent edge {:?} must exist for merge", info.parent_edge_key))?;
    let child_edge = self
      .edges
      .remove(&info.child_edge_key)
      .ok_or_else(|| make_internal_report!("Child edge {:?} must exist for merge", info.child_edge_key))?;

    self.nodes.remove(&info.removed_node_key);

    let merged_subs = parent_edge.chain_fitch_subs(child_edge.fitch_subs())?;
    let merged_indels = parent_edge.chain_fitch_indels(&child_edge.indels);

    let mut merged_edge = SparseEdgePartition::default();
    merged_edge.set_fitch_subs(merged_subs);
    merged_edge.indels = merged_indels;

    self.edges.insert(info.merged_edge_key, merged_edge);
    Ok(())
  }

  fn derive_root_sequence(&mut self, changes: &RerootChanges) {
    if !changes.inverted_edge_keys.is_empty() {
      let mut new_root_seq = self.root_sequence.clone();
      for edge_key in &changes.inverted_edge_keys {
        if let Some(edge_data) = self.edges.get(edge_key) {
          Self::apply_edge_to_sequence(&mut new_root_seq, edge_data, &self.alphabet);
        }
      }
      self.root_sequence = new_root_seq;
    } else if let Some(info) = &changes.edge_merge {
      if let Some(parent_edge) = self.edges.get(&info.parent_edge_key) {
        let parent_edge = parent_edge.clone();
        Self::apply_edge_to_sequence(&mut self.root_sequence, &parent_edge, &self.alphabet);
      }
    }
  }

  fn apply_edge_to_sequence(seq: &mut Seq, edge: &SparseEdgePartition, alphabet: &Alphabet) {
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
}

impl PartitionRerootOps for PartitionMarginalSparse {
  fn apply_reroot(&mut self, changes: &RerootChanges) -> Result<(), Report> {
    self.apply_reroot_changes(changes)?;

    if let Some(info) = &changes.edge_merge {
      self.remove_trivial_root(info)?;
    }

    Ok(())
  }
}
