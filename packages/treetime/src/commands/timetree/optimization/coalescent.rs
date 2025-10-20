use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::{GraphTimetree, PartitionTreetimeMarginalOps};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

/// Add coalescent model as population genetic prior on node times.
///
/// Optional: Improves time estimates by incorporating population dynamics.
/// Why: Coalescent theory provides realistic expectations for divergence times.
/// How: Kingman coalescent with exponential waiting times between mergers.
pub fn add_coalescent_model(
  _graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>],
  _params: &str,
) -> Result<(), Report> {
  todo!("Add merger rate contribution to node time distributions")
}
