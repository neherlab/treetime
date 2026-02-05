use crate::commands::clock::clock_model::ClockModel;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::GraphTimetree;
use eyre::Report;
use parking_lot::RwLock;
use std::path::Path;
use std::sync::Arc;

/// Assess sensitivity of timetree to clock rate uncertainty.
///
/// Optional: Quantifies robustness of time estimates.
/// Why: Clock rate has confidence intervals; propagating uncertainty improves reliability.
/// How: Re-run timetree with rate ± std_dev, compare resulting node times.
pub fn compute_rate_susceptibility(
  _graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  _clock_model: &ClockModel,
) -> Result<(), Report> {
  todo!("Run timetree with rate±σ, store alternative time estimates")
}

/// Extract confidence intervals from marginal posterior distributions.
///
/// Core: Quantifies uncertainty in inferred node times.
/// Why: Point estimates alone don't convey reliability of time inference.
/// How: Compute HPD (highest posterior density) or quantiles from marginal distributions.
pub fn extract_confidence_intervals(
  _graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
) -> Result<(), Report> {
  todo!("For each node: compute 95% HPD interval from marginal_pos_LH distribution")
}

/// Write confidence intervals for node dates.
///
/// Core: Essential for assessing reliability of time estimates.
/// Why: Quantifies uncertainty in molecular dating.
/// How: TSV file with node_name, lower_bound, median, upper_bound.
pub fn write_confidence_intervals(
  _graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTimetreeAll<NodeTimetree, EdgeTimetree>>>],
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Write confidence intervals to TSV file")
}
