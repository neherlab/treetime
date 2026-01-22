use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::GraphTimetree;
use eyre::Report;
use parking_lot::RwLock;
use std::path::Path;
use std::sync::Arc;

/// Write inferred node dates to file.
///
/// Core: Primary output of timetree inference.
/// Why: Users need calendar dates for each node.
/// How: TSV file with node_name, numdate, time_before_present.
pub fn write_node_dates(
  _graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Write node dates to TSV file")
}
