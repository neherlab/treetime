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
/// - **Dense**: returns MAP (maximum a posteriori) estimates from per-edge
///   probability profiles. MAP = argmax of each position's probability vector,
///   i.e. the single most probable state given the data and model.
pub trait PartitionOptimizeOps: Send + Sync {
  /// Return the sequence length represented by this partition.
  fn sequence_length(&self) -> usize;

  /// Return the precomputed likelihood contribution for one edge.
  fn create_edge_contribution(&self, edge_key: GraphEdgeKey) -> Result<OptimizationContribution, Report>;

  /// Return nucleotide substitutions for one edge.
  ///
  /// For sparse partitions, reconstructs parent and child states from Fitch
  /// parsimony data and current variable-site assignments, then compares.
  /// For dense partitions, takes the MAP state (argmax of the probability
  /// profile) at each end of the edge and compares. Gap positions are excluded
  /// in both cases.
  fn edge_subs(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<Vec<Sub>, Report>;

  /// Return the number of alignment positions where both parent and child
  /// have canonical (non-gap, non-ambiguous) states for one edge.
  fn edge_effective_length(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<usize, Report>;

  /// Return the sum of per-position differences for initial branch length estimation.
  ///
  /// For sparse partitions (default): integer count of discrete substitutions, equivalent
  /// to `edge_subs().len()`. For dense partitions: soft Hamming distance, where each
  /// position contributes `1 - dot(pp, pc)` (the complement of the profile overlap).
  /// Uncertain positions contribute fractionally rather than as 0 or 1.
  fn edge_initial_differences(&self, graph: &GraphAncestral, edge_key: GraphEdgeKey) -> Result<f64, Report> {
    Ok(self.edge_subs(graph, edge_key)?.len() as f64)
  }
}

pub type PartitionOptimizeVec = Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>;
