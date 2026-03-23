use crate::commands::optimize::optimize_unified::OptimizationContribution;
use crate::representation::payload::ancestral::GraphAncestral;
use crate::seq::mutation::Sub;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;

/// Operations that `optimize` needs from a sequence partition.
///
/// Dense and sparse partitions derive substitutions differently:
/// - **Sparse**: returns reconstructed states from Fitch parsimony and current
///   variable-site assignments (discrete, deterministic).
/// - **Dense**: returns MAP (maximum a posteriori) estimates from endpoint
///   node posteriors. After marginal inference, each node's `profile.dis`
///   contains the full posterior distribution over states. `edge_subs()`
///   compares the argmax state at the parent and child nodes.
pub trait PartitionOptimizeOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return the precomputed likelihood contribution for one edge.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Return nucleotide substitutions for one edge.
  ///
  /// For sparse partitions, reconstructs parent and child states from Fitch
  /// parsimony data and current variable-site assignments, then compares.
  /// For dense partitions, takes the MAP state (argmax of the node posterior)
  /// at the parent and child endpoints. Non-canonical states and gap positions
  /// are excluded in both cases.
  fn edge_subs(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Return the number of alignment positions where both parent and child
  /// have canonical (non-gap, non-ambiguous) states for one edge.
  fn edge_effective_length(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<usize, Report>;
}

pub type PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>;
