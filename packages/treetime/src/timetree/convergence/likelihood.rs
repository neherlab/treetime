use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::total_lh::compute_coalescent_total_lh;
use crate::progress::LogSink;
use crate::progress_warn;
use crate::timetree::branch_model::BranchModel;
use log::debug;
use std::collections::BTreeMap;
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::LogLh;

pub(crate) fn compute_sequence_log_lh(graph: &Graph, branch_model: &BranchModel) -> Option<LogLh> {
  let BranchModel::Marginal(partitions) = branch_model else {
    return None;
  };
  match partitions.graph_log_lh(graph) {
    Ok(lh) => Some(lh),
    Err(e) => {
      debug!("Sequence log-likelihood unavailable: {e}");
      None
    },
  }
}

pub(crate) fn compute_coalescent_log_lh(
  graph: &Graph,
  tc: &Distribution,
  node_times: &CoalescentNodeTimes,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  log: &dyn LogSink,
) -> Option<LogLh> {
  match compute_coalescent_total_lh(graph, tc, node_times, names, log) {
    Ok(lh) => Some(lh),
    Err(e) => {
      progress_warn!(log, "Coalescent log-likelihood unavailable: {e}");
      None
    },
  }
}
