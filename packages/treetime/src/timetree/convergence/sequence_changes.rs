use crate::seq::overlay::SeqOverlay;
use crate::timetree::branch_model::BranchModel;
use log::debug;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

pub(crate) fn count_sequence_changes(previous: &AncestralStateSnapshot, current: &AncestralStateSnapshot) -> usize {
  let prev_only = previous.keys().filter(|k| !current.contains_key(k)).count();
  let curr_only = current.keys().filter(|k| !previous.contains_key(k)).count();
  if prev_only > 0 || curr_only > 0 {
    debug!("{prev_only} nodes removed, {curr_only} nodes added between snapshots");
  }

  previous
    .iter()
    .filter_map(|(key, prev_seq)| current.get(key).map(|curr_seq| prev_seq.count_differences(curr_seq)))
    .sum()
}

pub(crate) fn capture_ancestral_states(graph: &Graph, branch_model: &BranchModel) -> AncestralStateSnapshot {
  match branch_model {
    BranchModel::Input => AncestralStateSnapshot::new(),
    BranchModel::Marginal(partition) => graph
      .get_nodes()
      .filter(|node| !node.is_leaf())
      .map(|node| (node.key(), partition.extract_ancestral_sequence(node.key())))
      .collect(),
  }
}

pub(crate) type AncestralStateSnapshot = BTreeMap<GraphNodeKey, SeqOverlay>;
