use crate::coalescent::events::collect_tree_events;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::coalescent::node_time::CoalescentNodeTimes;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

/// The number of extant lineages $k(t)$ implied by a tree's topology and node times. v0's
/// `nlineages`.
///
/// This is the fixed input of a coalescent model, as opposed to $T_c$, which is estimated, and
/// $H(t)$, which is a compound of the two. Which tree it is read from matters: $k(t)$ has two
/// roles, and only one of them may track the times being inferred. See
/// [`CoalescentModel`](crate::coalescent::coalescent::CoalescentModel).
pub fn compute_lineage_counts(graph: &Graph, node_times: &CoalescentNodeTimes) -> Result<PiecewiseConstantFn, Report> {
  let (_, events, terminal_lineage_count) = collect_tree_events(graph, node_times)?;
  compute_lineage_count_distribution(&events, terminal_lineage_count)
}
