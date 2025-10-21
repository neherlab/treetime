use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::graph::edge::GraphEdgeKey;
use crate::representation::graph_ancestral::{EdgeAncestral, NodeAncestral};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_marginal::PartitionMarginalOps;
use eyre::Report;

/// Trait for timetree-specific partition operations
/// Separate from PartitionMarginalOps to avoid polluting the ancestral command
pub trait PartitionTimetreeOps: Send + Sync {
  /// Create optimization contribution for branch length likelihood computation
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;
}

/// Combined trait for partitions that support both marginal and timetree operations
/// This allows trait objects to be used for both ancestral reconstruction and timetree inference
pub trait PartitionTimetreeAll:
  PartitionMarginalOps<NodeAncestral, EdgeAncestral> + PartitionTimetreeOps + HasLogLh
{
}

/// Blanket implementation: any type implementing all three traits automatically implements the combined trait
impl<T> PartitionTimetreeAll for T where
  T: PartitionMarginalOps<NodeAncestral, EdgeAncestral> + PartitionTimetreeOps + HasLogLh
{
}
