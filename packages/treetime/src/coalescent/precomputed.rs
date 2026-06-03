use crate::coalescent::events::collect_tree_events;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::coalescent::time_coordinate::CalendarTime;
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

/// Precomputed coalescent data derived from tree topology and node times.
///
/// Bundles the event collection and lineage count computation that every
/// coalescent entry point needs. Construct once via [`from_graph`](Self::from_graph),
/// then borrow the fields in downstream computations.
pub struct CoalescentPrecomputed {
  pub present_time: CalendarTime,
  pub lineage_counts: PiecewiseConstantFn,
}

impl CoalescentPrecomputed {
  pub fn from_graph<N, E, D>(graph: &Graph<N, E, D>) -> Result<Self, Report>
  where
    N: GraphNode + TimetreeNode,
    E: GraphEdge,
    D: Sync + Send,
  {
    let (present_time, events_calendar) = collect_tree_events(graph)?;
    let events_tbp: Vec<_> = events_calendar
      .iter()
      .map(|(t, delta)| (t.to_tbp(present_time), *delta))
      .collect();
    let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
    Ok(Self {
      present_time,
      lineage_counts,
    })
  }
}
