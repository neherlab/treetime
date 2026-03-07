use crate::commands::clock::reroot::RerootChanges;
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::partition::traits::PartitionMarginalOps;
use eyre::Report;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};

/// Trait for partition updates during reroot operations.
///
/// Default implementation is a no-op, suitable for dense partitions that don't track mutations.
/// Sparse partitions override to update mutation assignments when edges are split, merged, or inverted.
pub trait PartitionRerootOps: Send + Sync {
  /// Apply all partition updates for a reroot operation.
  ///
  /// Called after graph topology changes are complete. The `changes` struct bundles:
  /// - Edge split info (if a new node was created)
  /// - Edge merge info (if the old root was removed)
  /// - Inverted edge keys (path from old to new root)
  fn apply_reroot(&mut self, _changes: &RerootChanges) -> Result<(), Report> {
    Ok(())
  }
}

/// Trait for timetree-specific partition operations.
///
/// Separate from PartitionMarginalOps to avoid polluting the ancestral command.
pub trait PartitionTimetreeOps<N, E>: Send + Sync
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  /// Create optimization contribution for branch length likelihood computation.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Ensure partition has entries for all nodes and edges in the graph.
  ///
  /// After topology changes (polytomy resolution), new nodes/edges may lack partition entries.
  /// This adds empty/default entries for missing elements and removes stale entries for
  /// elements no longer in the graph. The subsequent marginal update pass recomputes values.
  fn reconcile_topology(&mut self, graph: &Graph<N, E, ()>);
}

/// Combined trait for partitions that support both marginal and timetree operations.
/// This allows trait objects to be used for both ancestral reconstruction and timetree inference.
/// Includes `PartitionRerootOps` for reroot support with default no-op.
pub trait PartitionTimetreeAll<N, E>:
  PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + PartitionRerootOps + HasLogLh
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
}

/// Blanket implementation: any type implementing all required traits automatically implements the combined trait
impl<T, N, E> PartitionTimetreeAll<N, E> for T
where
  T: PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + PartitionRerootOps + HasLogLh,
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
}
