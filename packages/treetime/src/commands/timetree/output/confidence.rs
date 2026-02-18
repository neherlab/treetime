use crate::commands::clock::clock_model::ClockModel;
use crate::commands::timetree::partition_ops::PartitionTimetreeAll;
use crate::representation::partition::timetree::GraphTimetree;
use crate::representation::payload::timetree::EdgeTimetree;
use crate::representation::payload::timetree::NodeTimetree;
use eyre::Report;
use itertools::Itertools;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::path::Path;
use std::sync::Arc;
use treetime_graph::node::{Named, TimeConstraint};
use treetime_io::csv::CsvStructFileWriter;

/// Confidence interval for a node's inferred date.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NodeConfidenceInterval {
  pub name: String,
  pub date: f64,
  pub lower: f64,
  pub upper: f64,
}

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
/// Computes 95% CI for each node from time distribution using quantiles at 0.05 and 0.95.
/// For delta distributions (point estimates), returns identity interval `[date, date]`.
pub fn extract_confidence_intervals(graph: &GraphTimetree) -> Vec<NodeConfidenceInterval> {
  graph
    .get_nodes()
    .into_iter()
    .filter_map(|node_ref| {
      let node = node_ref.read_arc();
      let payload = node.payload().read_arc();
      let name = payload.name().map(|n| n.as_ref().to_owned())?;
      let date = payload.time?;

      let (lower, upper) = payload
        .time_distribution()
        .as_ref()
        .and_then(|dist| dist.confidence_interval(0.05, 0.95))
        .unwrap_or((date, date));

      Some(NodeConfidenceInterval {
        name,
        date,
        lower,
        upper,
      })
    })
    .sorted_by(|a, b| a.name.cmp(&b.name))
    .collect_vec()
}

/// Write confidence intervals for node dates to TSV file.
pub fn write_confidence_intervals(intervals: &[NodeConfidenceInterval], filepath: &Path) -> Result<(), Report> {
  let mut writer = CsvStructFileWriter::new(filepath, b'\t')?;
  intervals.iter().try_for_each(|ci| writer.write(ci))
}

/// Combine two confidence interval contributions via quadrature sum.
///
/// Given a center date and physical limits, combines up to two independent
/// confidence interval contributions (e.g., from mutation uncertainty and rate uncertainty).
/// Returns bounds clipped to the physical limits.
///
/// Formula:
/// - `lower = center - sqrt((c1_lo - center)² + (c2_lo - center)²)`
/// - `upper = center + sqrt((c1_hi - center)² + (c2_hi - center)²)`
pub fn combine_confidence(
  center: f64,
  limits: (f64, f64),
  c1: Option<(f64, f64)>,
  c2: Option<(f64, f64)>,
) -> (f64, f64) {
  let (min_val, max_val) = match (c1, c2) {
    (None, None) => return limits,
    (Some(c), None) | (None, Some(c)) => c,
    (Some(c1), Some(c2)) => {
      let min_val = center - (c1.0 - center).hypot(c2.0 - center);
      let max_val = center + (c1.1 - center).hypot(c2.1 - center);
      (min_val, max_val)
    },
  };

  (limits.0.max(min_val), limits.1.min(max_val))
}
