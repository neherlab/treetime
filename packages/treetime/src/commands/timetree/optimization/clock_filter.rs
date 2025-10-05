use crate::commands::timetree::data::date_constraints::DateConstraintSet;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;

/// Detect and mark branches that violate molecular clock assumptions.
///
/// Optional: Improves robustness by excluding outliers from clock model inference.
/// Why: Recombination, selection, or sequencing errors can create non-clock-like branches.
/// How: Root-to-tip regression residuals analyzed via IQD (interquartile distance) threshold.
pub fn clock_filter(_graph: &GraphAncestral, _constraints: &DateConstraintSet, _n_iqd: f64) -> Result<(), Report> {
  todo!("Mark branches as bad_branch=true if |residual| > n_iqd * IQD")
}

/// Identify and report outlier branches that violate molecular clock.
///
/// Core: User needs to know which samples were excluded from analysis.
/// Why: Outliers may indicate data quality issues requiring investigation.
/// How: Compare inferred dates to input constraints, report significant deviations.
pub fn report_bad_branches(_graph: &GraphAncestral, _constraints: &DateConstraintSet) -> Result<(), Report> {
  todo!("Log warnings for branches where |inferred_date - constraint_date| is large")
}
