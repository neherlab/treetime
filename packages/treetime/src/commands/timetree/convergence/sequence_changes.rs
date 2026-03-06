use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

/// Per-partition map of internal node keys to their reconstructed sequences.
pub type AncestralStateSnapshot = Vec<BTreeMap<GraphNodeKey, Seq>>;

/// Count positions where ancestral sequences differ between two snapshots.
///
/// Compares position-by-position across all internal nodes and partitions.
pub fn count_sequence_changes(previous: &AncestralStateSnapshot, current: &AncestralStateSnapshot) -> usize {
  previous
    .iter()
    .zip(current.iter())
    .map(|(prev_partition, curr_partition)| {
      prev_partition
        .iter()
        .filter_map(|(key, prev_seq)| {
          curr_partition
            .get(key)
            .map(|curr_seq| count_differing_positions(prev_seq, curr_seq))
        })
        .sum::<usize>()
    })
    .sum()
}

/// Snapshot current ancestral sequences for all internal nodes across partitions.
pub fn capture_ancestral_states(
  graph: &GraphTimetree,
  partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
) -> AncestralStateSnapshot {
  if partitions.is_empty() {
    return vec![];
  }

  let internal_keys: Vec<GraphNodeKey> = graph.filter_map(|node| if node.is_leaf { None } else { Some(node.key) });

  partitions
    .iter()
    .map(|partition| {
      let partition = partition.read_arc();
      internal_keys
        .iter()
        .map(|&key| (key, partition.extract_ancestral_sequence(key)))
        .collect()
    })
    .collect()
}

fn count_differing_positions(a: &Seq, b: &Seq) -> usize {
  let shared = a.iter().zip(b.iter()).filter(|(ca, cb)| ca != cb).count();
  let length_diff = a.len().abs_diff(b.len());
  shared + length_diff
}
