use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::representation::partition::traits::PartitionBranchOps;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;

/// Operations that `optimize` needs from a sequence partition.
///
/// Extends `PartitionBranchOps` (which provides `edge_subs()`,
/// `edge_effective_length()`, and `sequence_length()`) with
/// optimize-specific likelihood contribution computation.
pub trait PartitionOptimizeOps: PartitionBranchOps + Send + Sync {
  /// Return the precomputed likelihood contribution for one edge.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;
}

pub type PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>;
