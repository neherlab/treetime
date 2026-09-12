use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::total_lh::compute_coalescent_total_lh;
use crate::partition::timetree::partition::PartitionTimetree;
use crate::partition::timetree::marginal::graph_log_lh;
use crate::timetree::timetree_state::TimetreeState;
use log::{debug, warn};
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;

/// Sum of per-partition root log-likelihoods from marginal reconstruction.
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

/// Sum of log-probabilities of branch length distributions evaluated at inferred time durations.
///
/// For each edge, evaluates the branch length distribution at `child_time - parent_time`
/// (positive duration in calendar time where parent < child).
///
/// This is a v1-specific metric. v0's `positional_LH` sums node-level marginal
/// log-likelihoods from the forward pass. Both metrics trend in the same direction
/// during convergence but produce different numerical values.
pub fn compute_positional_log_lh(graph: &Graph, state: &TimetreeState) -> Option<LogLh> {
  let mut total = 0.0;
  let mut count = 0_usize;

  for edge_ref in graph.get_edges() {
    let edge = edge_ref.read_arc();
    let parent_key = edge.source();
    let child_key = edge.target();

    let edge_state = state.edge(edge.key());
    let Some(dist) = edge_state.branch_length_distribution.as_ref() else {
      continue;
    };

    let parent_time = state.node(parent_key).time;
    let child_time = state.node(child_key).time;

    // Calendar time: parent_time < child_time, so duration = ct - pt is positive,
    // matching the positive domain of branch length distributions.
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

/// Total coalescent log-likelihood from the merger model.
///
/// Sums per-edge costs under the Kingman coalescent for the given Tc distribution.
/// Returns `None` when no coalescent model is active (coalescent_tc is None).
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
