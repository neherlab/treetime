use crate::alphabet::alphabet::Alphabet;
use crate::partition::assembly::AlignmentGroup;
use crate::partition::marginal::reconstruction::SparseReconstruction;
use crate::partition::marginal::sparse::alignment_indels::{
  assembled_edge_indels, merge_edge_pair, sparse_alignment_groups, store_edge_indels,
};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
use crate::{make_internal_error, make_internal_report};
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::{EdgeMergeInfo, RerootResult, StemRemovalInfo};
use treetime_primitives::Seq;

pub(crate) fn reroot_sparse(
  partitions: Vec<SparseReconstruction>,
  changes: &RerootResult,
) -> Result<Vec<SparseReconstruction>, Report> {
  let groups = sparse_alignment_groups(&partitions)?;
  let mut partitions = partitions;
  for SparseReconstruction {
    partition, node_states, ..
  } in &mut partitions
  {
    if let Some(info) = &changes.stem_removal {
      remove_stem_root(partition, node_states, info)?;
    }
    split_reroot_edge(partition, changes)?;
  }
  for edge_key in &changes.inverted_edge_keys {
    invert_edge(&mut partitions, &groups, *edge_key)?;
  }
  for SparseReconstruction {
    partition, node_states, ..
  } in &mut partitions
  {
    derive_root_sequence(partition, changes);
    if let Some(info) = &changes.edge_split {
      partition.obs_nodes.insert(
        info.new_node_key,
        SparseNodeObs::new(&partition.root_sequence, &partition.alphabet),
      );
      node_states.insert(info.new_node_key, SparseNodeState::leaf(&partition.root_sequence));
    }
  }
  if let Some(info) = &changes.edge_merge {
    remove_trivial_root(&mut partitions, &groups, info)?;
  }
  Ok(
    partitions
      .into_iter()
      .map(
        |SparseReconstruction {
           partition,
           gtr,
           node_states,
           layout,
           ..
         }| {
          let node_states = reconcile_node_states(&partition, node_states);
          SparseReconstruction::seeded(partition, gtr, node_states).keeping_layout(layout)
        },
      )
      .collect(),
  )
}

fn split_reroot_edge(partition: &mut PartitionMarginalSparse, changes: &RerootResult) -> Result<(), Report> {
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
  Ok(())
}

fn invert_edge(
  partitions: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  edge_key: GraphEdgeKey,
) -> Result<(), Report> {
  for SparseReconstruction { partition, .. } in partitions.iter_mut() {
    partition
      .obs_edges
      .get_mut(&edge_key)
      .ok_or_else(|| make_internal_report!("Edge {edge_key:?} must exist on reroot path"))?
      .invert_fitch_subs();
  }
  for group in groups {
    let mut indels = assembled_edge_indels(group, partitions, edge_key)?;
    for indel in &mut indels {
      indel.invert();
    }
    store_edge_indels(group, partitions, edge_key, &indels)?;
  }
  Ok(())
}

fn remove_stem_root(
  partition: &mut PartitionMarginalSparse,
  node_states: &mut BTreeMap<GraphNodeKey, SparseNodeState>,
  info: &StemRemovalInfo,
) -> Result<(), Report> {
  let mut stem_edge = partition
    .obs_edges
    .remove(&info.removed_edge_key)
    .ok_or_else(|| make_internal_report!("Stem edge {:?} must exist for removal", info.removed_edge_key))?;
  stem_edge.invert_fitch_subs();
  for indel in &mut stem_edge.indels {
    indel.invert();
  }
  apply_edge_to_sequence(&mut partition.root_sequence, &stem_edge, &partition.alphabet);
  partition.obs_nodes.remove(&info.removed_node_key);
  node_states.remove(&info.removed_node_key);
  Ok(())
}

fn remove_trivial_root(
  partitions: &mut [SparseReconstruction],
  groups: &[AlignmentGroup],
  info: &EdgeMergeInfo,
) -> Result<(), Report> {
  for SparseReconstruction {
    partition, node_states, ..
  } in partitions.iter_mut()
  {
    for (edge_key, role) in [(info.parent_edge_key, "Parent"), (info.child_edge_key, "Child")] {
      if !partition.obs_edges.contains_key(&edge_key) {
        return make_internal_error!("{role} edge {edge_key:?} must exist for merge");
      }
    }
    partition.obs_nodes.remove(&info.removed_node_key);
    node_states.remove(&info.removed_node_key);
  }
  merge_edge_pair(partitions, groups, info)
}

fn derive_root_sequence(partition: &mut PartitionMarginalSparse, changes: &RerootResult) {
  for edge_key in changes.inverted_edge_keys.iter().rev() {
    apply_edge_to_sequence(
      &mut partition.root_sequence,
      &partition.obs_edges[edge_key],
      &partition.alphabet,
    );
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
