use crate::representation::partition_timetree::GraphTimetree;
use eyre::Report;

/// Detect and mark branches that violate molecular clock assumptions.
///
/// Optional: Improves robustness by excluding outliers from clock model inference.
/// Why: Recombination, selection, or sequencing errors can create non-clock-like branches.
/// How: Root-to-tip regression residuals analyzed via IQD (interquartile distance) threshold.
pub fn filter_outliers(_graph: &GraphTimetree, _n_iqd: f64) -> Result<(), Report> {
  todo!("filter_outliers")
}

pub fn mark_bad_branches(_graph: &GraphTimetree, _n_iqd: f64) -> Result<(), Report> {
  todo!("mark_bad_branches")
}

pub fn clock_filter(_graph: &GraphTimetree, _n_iqd: f64) -> Result<(), Report> {
  todo!("clock_filter")
}

pub fn reroot_on_best(_graph: &GraphTimetree, _resolve_polytomies: bool) -> Result<(), Report> {
  todo!("reroot_on_best")
}

/// Identify and report outlier branches that violate molecular clock.
///
/// Core: User needs to know which samples were excluded from analysis.
/// Why: Outliers may indicate data quality issues requiring investigation.
/// How: Compare inferred dates to input constraints, report significant deviations.
pub fn report_bad_branches(_graph: &GraphTimetree) -> Result<(), Report> {
  todo!("report_bad_branches")
}
