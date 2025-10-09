use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_timetree::PartitionTreetimeMarginalOps;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::sync::Arc;

/// Stores ancestral sequence states for change tracking between iterations.
///
/// Maps node keys to their reconstructed sequences to enable ndiff calculation.
pub type AncestralStateSnapshot = BTreeMap<String, Vec<u8>>;

/// Count ancestral sequence changes between iterations.
///
/// Core: Primary convergence criterion (stop when ndiff == 0).
/// Why: No sequence changes indicates reconstruction has stabilized.
/// How: Compare current and previous ancestral sequences, count differing positions.
pub fn count_sequence_changes(
  _previous_states: &AncestralStateSnapshot,
  _current_states: &AncestralStateSnapshot,
) -> usize {
  // TODO: Compare sequences position-by-position across all internal nodes
  // Return total count of positions where ancestral state changed
  0
}

/// Capture current ancestral sequence states for change tracking.
///
/// Core: Required for ndiff calculation between iterations.
/// Why: Need snapshot to compare before/after reconstruction.
/// How: Extract reconstructed sequences from all internal nodes.
pub fn capture_ancestral_states(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
) -> AncestralStateSnapshot {
  // TODO: Extract ancestral sequences from partition node data
  // Map node keys to their current sequence states
  BTreeMap::new()
}
