use crate::coalescent::events::collect_tree_events;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::coalescent::node_time::CoalescentNodeTimes;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

pub fn compute_lineage_counts(graph: &Graph, node_times: &CoalescentNodeTimes) -> Result<PiecewiseConstantFn, Report> {
  let (_, events, terminal_lineage_count) = collect_tree_events(graph, node_times)?;
  compute_lineage_count_distribution(&events, terminal_lineage_count)
}
