use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_timetree::PartitionTreetimeMarginalOps;
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

/// Add coalescent model as population genetic prior on node times.
///
/// Optional: Improves time estimates by incorporating population dynamics.
/// Why: Coalescent theory provides realistic expectations for divergence times.
/// How: Kingman coalescent with exponential waiting times between mergers.
pub fn add_coalescent_model(
  _graph: &GraphAncestral,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps>>],
  _params: &str,
) -> Result<(), Report> {
  todo!("Add merger rate contribution to node time distributions")
}
