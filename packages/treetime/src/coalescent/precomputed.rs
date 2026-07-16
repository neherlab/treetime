use crate::coalescent::events::collect_tree_events;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
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
#[derive(Clone, Debug)]
pub struct CoalescentPrecomputed {
  lineage_counts: PiecewiseConstantFn,
}

impl CoalescentPrecomputed {
  pub fn lineage_counts(&self) -> &PiecewiseConstantFn {
    &self.lineage_counts
  }

  #[cfg(test)]
  pub fn from_lineage_counts(lineage_counts: PiecewiseConstantFn) -> Self {
    Self { lineage_counts }
  }

  pub fn from_graph<N, E, D>(graph: &Graph<N, E, D>) -> Result<Self, Report>
  where
    N: GraphNode + TimetreeNode,
    E: GraphEdge,
    D: Sync + Send,
  {
    let (_, events) = collect_tree_events(graph)?;
    let lineage_counts = compute_lineage_count_distribution(&events)?;
    Ok(Self { lineage_counts })
  }
}
