use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::total_lh::compute_coalescent_total_lh;
use crate::progress::LogSink;
use crate::progress_warn;
use crate::timetree::branch_model::BranchModel;
use crate::timetree::inference::time_inference::TimeInference;
use log::debug;
use std::collections::BTreeMap;
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

pub(crate) fn compute_sequence_log_lh(graph: &Graph, branch_model: &BranchModel) -> Option<LogLh> {
  let BranchModel::Marginal(partition) = branch_model else {
    return None;
  };
  match partition.graph_log_lh(graph) {
    Ok(lh) => Some(lh),
    Err(e) => {
      debug!("Sequence log-likelihood unavailable: {e}");
      None
    },
  }
}

pub(crate) fn compute_positional_log_lh(graph: &Graph, inference: &TimeInference) -> Option<LogLh> {
  let mut total = 0.0;
  let mut count = 0_usize;

  for edge_ref in graph.get_edges() {
    let edge = edge_ref;
    let parent_key = edge.source();
    let child_key = edge.target();

    let Some(dist) = inference.branches[&edge.key()].distribution.as_ref() else {
      continue;
    };

    let parent_time = inference.posterior[&parent_key].time;
    let child_time = inference.posterior[&child_key].time;

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

pub(crate) fn compute_coalescent_log_lh(
  graph: &Graph,
  coalescent_tc: Option<&Distribution>,
  node_times: &CoalescentNodeTimes,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Option<LogLh> {
  let tc = coalescent_tc?;
  match compute_coalescent_total_lh(graph, tc, node_times, names, log) {
    Ok(lh) => Some(lh),
    Err(e) => {
      progress_warn!(log, "Coalescent log-likelihood unavailable: {e}");
      None
    },
  }
}
