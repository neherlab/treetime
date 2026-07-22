use crate::partition::timetree::{GraphTimetree, PartitionTimetreeRef};
use eyre::Report;
use std::path::Path;

/// Write inferred node dates to file.
///
/// Core: Primary output of timetree inference.
/// Why: Users need calendar dates for each node.
/// How: TSV file with node_name, numdate, time_before_present.
pub fn write_node_dates(
  _graph: &GraphTimetree,
  _partitions: &[PartitionTimetreeRef],
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Write node dates to TSV file")
}
