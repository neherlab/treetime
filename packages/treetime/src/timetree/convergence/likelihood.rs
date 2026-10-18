use crate::timetree::branch_model::BranchModel;
use log::debug;
use treetime_graph::graph::Graph;
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
