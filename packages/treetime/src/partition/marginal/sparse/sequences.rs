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
  let inheritance = &partition.inheritance;
  let nodes = node_states.iter().map(|(key, node)| (key, &node.sequence)).collect();
  let edges = partition.obs_edges.iter().collect();
  let outputs = schedule.scan_down(
    inheritance.join_parent_edge(schedule, nodes, edges),
    |key, (sequence, edge_obs), parents: PerEdge<Option<Arc<Seq>>>| {
      let Some(sequence) = sequence else {
        return Ok((None, None));
      };
      let keeps_sequence = inheritance.children(schedule, key).next().is_none() && !sequence.is_empty();
      let parent = inheritance.parent_message(schedule, key, &parents)?;
      let sequence = match (edge_obs, parent) {
        (Some(_), Some(_)) if keeps_sequence => Arc::clone(sequence),
        (Some(edge_obs), Some((_, _, parent_sequence))) => {
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
      Ok((Some(Arc::clone(&sequence)), Some(sequence)))
    },
  )?;
  Ok(
    outputs
      .into_map()
      .into_iter()
      .filter_map(|(key, sequence)| sequence.map(|sequence| (key, sequence)))
      .collect(),
  )
}
