use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::marginal::sparse::reconstruct::parsimony_seq;
use crate::partition::storage::sparse::SparseNodeState;
use eyre::Report;
use std::sync::Arc;
use treetime_graph::key_map::NodeMap;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::PerEdge;
use treetime_primitives::Seq;
use treetime_utils::make_internal_report;

pub(crate) fn fitch_node_sequences(
  partition: &PartitionMarginalSparse,
  schedule: &GraphSchedule<'_>,
  node_states: &NodeMap<SparseNodeState>,
) -> Result<NodeMap<Arc<Seq>>, Report> {
  let nodes = node_states.iter().map(|(key, node)| (key, &node.sequence)).collect();
  let edges = partition.obs_edges.iter().collect();
  let outputs = schedule.scan_down(
    schedule.join_parent_edges(nodes, edges)?,
    |key, (sequence, mut edge_obs), parents: PerEdge<Arc<Seq>>| {
      let keeps_sequence = schedule.children(key).is_empty() && !sequence.is_empty();
      let sequence = match (edge_obs.pop(), parents.first()) {
        (Some(_), Some(_)) if keeps_sequence => Arc::clone(sequence),
        (Some(edge_obs), Some(parent_sequence)) => {
          let rebuilt = parsimony_seq(
            parent_sequence,
            edge_obs,
            &partition.obs_nodes[&key],
            &partition.alphabet,
          );
          if **sequence == rebuilt {
            Arc::clone(sequence)
          } else {
            Arc::new(rebuilt)
          }
        },
        (None, None) if sequence.is_empty() => Arc::new(partition.root_sequence.clone()),
        (None, None) => Arc::clone(sequence),
        _ => {
          return Err(make_internal_report!(
            "Sparse sequence pass received an inconsistent parent edge and parent sequence at node {key}"
          ));
        },
      };
      Ok((Arc::clone(&sequence), sequence))
    },
  )?;
  Ok(outputs.into_map())
}
