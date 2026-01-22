use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::graph::edge::{GraphEdge, GraphEdgeKey, Weighted};
use crate::graph::node::{GraphNode, Named};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_marginal::PartitionMarginalOps;
use eyre::Report;

/// Trait for timetree-specific partition operations
/// Separate from PartitionMarginalOps to avoid polluting the ancestral command
pub trait PartitionTimetreeOps<N, E>: Send + Sync
where
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
{
  /// Create optimization contribution for branch length likelihood computation
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;
}

/// Combined trait for partitions that support both marginal and timetree operations
/// This allows trait objects to be used for both ancestral reconstruction and timetree inference
pub trait PartitionTimetreeAll<N, E>: PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + HasLogLh
where
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
{
}

/// Blanket implementation: any type implementing all three traits automatically implements the combined trait
impl<T, N, E> PartitionTimetreeAll<N, E> for T
where
  T: PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + HasLogLh,
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
{
}
