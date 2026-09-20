use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::total_lh::compute_coalescent_total_lh;
use crate::partition::timetree::marginal::graph_log_lh;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::timetree::timetree_state::TimetreeState;
use log::{debug, warn};
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;

pub fn compute_sequence_log_lh(graph: &Graph, partitions: &[PartitionTimetree]) -> Option<LogLh> {
  if partitions.is_empty() {
    return None;
  }
  match graph_log_lh(graph, partitions) {
    Ok(lh) => Some(lh),
    Err(e) => {
      debug!("Sequence log-likelihood unavailable: {e}");
      None
    },
  }
}

pub fn compute_positional_log_lh(graph: &Graph, state: &TimetreeState) -> Option<LogLh> {
  let mut total = 0.0;
  let mut count = 0_usize;

  for edge_ref in graph.get_edges() {
    let edge = edge_ref;
    let parent_key = edge.source();
    let child_key = edge.target();

    let edge_state = state.edge(edge.key());
    let Some(dist) = edge_state.branch_length_distribution.as_ref() else {
      continue;
    };

    let parent_time = state.node(parent_key).time;
    let child_time = state.node(child_key).time;

    let time_diff = match (parent_time, child_time) {
      (Some(pt), Some(ct)) => ct - pt,
      _ => continue,
    };

    match dist.eval(time_diff) {
      Ok(p) if p > 0.0 => {
        total += p.ln();
        count += 1;
      },
      Ok(p) => {
        debug!("Edge {parent_key:?}->{child_key:?}: zero or negative probability {p:.6e} at time_diff={time_diff:.6e}");
      },
      Err(e) => {
        debug!("Edge {parent_key:?}->{child_key:?}: distribution eval failed at time_diff={time_diff:.6e}: {e}");
      },
    }
  }

  (count > 0).then_some(LogLh::new(total))
}

pub fn compute_coalescent_log_lh(
  graph: &Graph,
  coalescent_tc: Option<&Distribution>,
  node_times: &CoalescentNodeTimes,
) -> Option<LogLh> {
  let tc = coalescent_tc?;
  match compute_coalescent_total_lh(graph, tc, node_times) {
    Ok(lh) => Some(lh),
    Err(e) => {
      warn!("Coalescent log-likelihood unavailable: {e}");
      None
    },
  }
}
