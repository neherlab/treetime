use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::representation::payload::ancestral::GraphAncestral;
use crate::seq::mutation::Sub;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;

/// Operations that `optimize` needs from a sequence partition.
///
/// In this trait, branch mutations mean the current nucleotide changes between
/// the parent node and the child node on one edge. Dense and sparse partitions
/// must return the same kind of result.
pub trait PartitionOptimizeOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return the precomputed likelihood contribution for one edge.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Return the current nucleotide changes for one edge.
  ///
  /// Sparse needs `graph` to find the edge endpoints and walk back to the
  /// parent state when needed. Dense ignores `graph` because it already has
  /// enough per-edge data.
  fn edge_subs(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;
}

pub type PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>;
