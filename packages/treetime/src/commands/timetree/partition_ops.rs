use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, GraphNodeKey, Named};
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_marginal::PartitionMarginalOps;
use eyre::Report;

/// Trait for timetree-specific partition operations
/// Separate from PartitionMarginalOps to avoid polluting the ancestral command
pub trait PartitionTimetreeOps<N, E>: Send + Sync
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  /// Create optimization contribution for branch length likelihood computation
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Whether this partition supports edge splitting during reroot
  fn supports_reroot_edge_split(&self) -> bool;

  /// Whether this partition supports old root removal during reroot
  fn supports_old_root_removal(&self) -> bool;

  /// Update partition state after a node-only reroot (no edge splitting)
  ///
  /// `path_from_new_to_old` is the path from new root to old root, matching
  /// the format returned by `Graph::path_from_node_to_node(new_root, old_root)`.
  /// Each element is `(node_key, Option<edge_key>)` where the edge connects
  /// to the next node in the path.
  fn reroot_partition_node_only(
    &mut self,
    graph: &Graph<N, E, ()>,
    old_root_key: GraphNodeKey,
    new_root_key: GraphNodeKey,
    path_from_new_to_old: &[(GraphNodeKey, Option<GraphEdgeKey>)],
  ) -> Result<(), Report>;
}

/// Combined trait for partitions that support both marginal and timetree operations
/// This allows trait objects to be used for both ancestral reconstruction and timetree inference
pub trait PartitionTimetreeAll<N, E>: PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + HasLogLh
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
}

/// Blanket implementation: any type implementing all three traits automatically implements the combined trait
impl<T, N, E> PartitionTimetreeAll<N, E> for T
where
  T: PartitionMarginalOps<N, E> + PartitionTimetreeOps<N, E> + HasLogLh,
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
}
