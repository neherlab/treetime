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
  /// Both sparse and dense use `graph` to find the edge endpoints. Dense
  /// also uses it to look up parent/child gap ranges.
  fn edge_subs(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Return the number of alignment positions where both parent and child
  /// have canonical (non-gap, non-ambiguous) states for one edge.
  fn edge_effective_length(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<usize, Report>;
}

pub type PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>;
