use crate::coalescent::coalescent::CoalescentModel;
use crate::coalescent::edge_data::{coalescent_log_likelihood, collect_coalescent_edges};
use crate::coalescent::lineage_counts::compute_lineage_counts;
use crate::coalescent::node_time::CoalescentNodeTimes;
use eyre::Report;
use treetime_distribution::Distribution;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;

pub(crate) fn compute_coalescent_total_lh(
  graph: &Graph,
  tc_dist: &Distribution,
  node_times: &CoalescentNodeTimes,
) -> Result<LogLh, Report> {
  let model = CoalescentModel::new(&compute_lineage_counts(graph, node_times)?, tc_dist)?;
  let edges = collect_coalescent_edges(graph, node_times)?;

  coalescent_log_likelihood(&edges, &model)
}
