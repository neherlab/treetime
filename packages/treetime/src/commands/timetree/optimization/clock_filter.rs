use crate::representation::partition::timetree::GraphTimetree;
use eyre::Report;

/// Identify and report outlier branches that violate molecular clock.
///
/// Core: User needs to know which samples were excluded from analysis.
/// Why: Outliers may indicate data quality issues requiring investigation.
/// How: Compare inferred dates to input constraints, report significant deviations.
pub fn report_bad_branches(_graph: &GraphTimetree) -> Result<(), Report> {
  todo!("report_bad_branches")
}
