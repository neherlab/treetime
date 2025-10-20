use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::PartitionTreetimeMarginalOps;
use parking_lot::RwLock;
use std::sync::Arc;

/// Calculate total sequence likelihood across all partitions.
///
/// Core: Component of convergence monitoring.
/// Why: Track changes in sequence reconstruction quality across iterations.
/// How: Sum sequence likelihood from all partitions (joint or marginal depending on mode).
pub fn calculate_sequence_likelihood(_partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>]) -> Option<f64> {
  // TODO: Extract sequence_joint_LH or sequence_marginal_LH from partitions
  // For now return None until partition likelihood tracking is implemented
  None
}

/// Calculate total positional (temporal) likelihood across all partitions.
///
/// Core: Component of convergence monitoring.
/// Why: Track changes in node time inference quality across iterations.
/// How: Sum positional likelihood from all partitions.
pub fn calculate_positional_likelihood(_partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>]) -> Option<f64> {
  // TODO: Extract positional_LH from partitions
  // This represents the likelihood of the inferred node times given temporal constraints
  None
}

/// Calculate total coalescent likelihood across all partitions.
///
/// Optional: Only relevant when coalescent model is active.
/// Why: Track coalescent prior contribution to convergence.
/// How: Sum coalescent likelihood from merger model if present.
pub fn calculate_coalescent_likelihood(_partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>]) -> Option<f64> {
  // TODO: Extract coalescent_joint_LH from merger model if active
  // Returns None when coalescent model not used
  None
}

/// Calculate total likelihood (sum of all components).
///
/// Core: Primary convergence metric.
/// Why: Combined likelihood should increase or stabilize as tree converges.
/// How: Sum sequence, positional, and coalescent likelihoods.
pub fn calculate_total_likelihood(partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>]) -> Option<f64> {
  let seq_lh = calculate_sequence_likelihood(partitions);
  let pos_lh = calculate_positional_likelihood(partitions);
  let coal_lh = calculate_coalescent_likelihood(partitions);

  match (seq_lh, pos_lh, coal_lh) {
    (Some(s), Some(p), Some(c)) => Some(s + p + c),
    (Some(s), Some(p), None) => Some(s + p),
    (Some(s), None, None) => Some(s),
    (None, Some(p), None) => Some(p),
    _ => None,
  }
}
