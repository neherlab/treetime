use crate::ancestral::pipeline::DenseReconstruction;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::storage::dense::DenseNodeState;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::RerootChanges;

pub(crate) fn reroot_dense(
  partition: PartitionMarginalDense,
  gtr: GTR,
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

  DenseReconstruction::seeded(partition, gtr, node_states)
}
