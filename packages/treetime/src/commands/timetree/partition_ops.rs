use crate::commands::clock::reroot::{EdgeMergeInfo, EdgeSplitInfo};
use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
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

  /// Handle edge split during reroot: create new node entry, move mutations to child-side edge
  fn handle_edge_split(&mut self, info: &EdgeSplitInfo) -> Result<(), Report>;

  /// Handle edge merge when removing trivial node: compose mutations from two edges
  fn handle_edge_merge(&mut self, info: &EdgeMergeInfo) -> Result<(), Report>;

  /// Update partition state after reroot path inversion
  ///
  /// `path_from_old_to_new` is the path from old root to new root (upward in post-reroot tree),
  /// matching the format returned by `Graph::path_from_node_to_node(old_root, new_root)`.
  /// Each element is `(node_key, Option<edge_key>)` where the edge connects
  /// to the next node in the path.
  fn update_partition_after_reroot(
    &mut self,
    old_root_key: GraphNodeKey,
    new_root_key: GraphNodeKey,
    path_from_old_to_new: &[(GraphNodeKey, Option<GraphEdgeKey>)],
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
