use crate::partition::timetree::partition::PartitionTimetree;
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::Seq;

pub type AncestralStateSnapshot = Vec<BTreeMap<GraphNodeKey, Seq>>;

pub fn count_sequence_changes(previous: &AncestralStateSnapshot, current: &AncestralStateSnapshot) -> usize {
  previous
    .iter()
    .zip(current.iter())
    .enumerate()
    .map(|(partition_idx, (prev_partition, curr_partition))| {
      let prev_only = prev_partition
        .keys()
        .filter(|k| !curr_partition.contains_key(k))
        .count();
      let curr_only = curr_partition
        .keys()
        .filter(|k| !prev_partition.contains_key(k))
        .count();
      if prev_only > 0 || curr_only > 0 {
        debug!("Partition {partition_idx}: {prev_only} nodes removed, {curr_only} nodes added between snapshots");
      }

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

pub fn capture_ancestral_states(graph: &Graph, partitions: &[PartitionTimetree]) -> AncestralStateSnapshot {
  if partitions.is_empty() {
    return vec![];
  }

  let internal_keys: Vec<GraphNodeKey> = graph
    .get_nodes()
    .filter(|node| !node.is_leaf())
    .map(|node| node.key())
    .collect();

  partitions
    .iter()
    .map(|partition| {
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
