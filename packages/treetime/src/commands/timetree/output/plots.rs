use crate::commands::timetree::data::date_constraints::DateConstraintSet;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;
use std::path::Path;

/// Plot root-to-tip distance vs sampling date regression.
///
/// Optional: Visual quality control for temporal signal.
/// Why: Diagnostic for clock-like evolution and outlier detection.
/// How: Scatter plot with regression line, residuals, R² annotation.
pub fn plot_root_to_tip(
  _graph: &GraphAncestral,
  _constraints: &DateConstraintSet,
  _out_base: &Path,
) -> Result<(), Report> {
  todo!("Generate root-to-tip regression plot")
}

/// Plot time-scaled phylogenetic tree.
///
/// Optional: Visual representation of inferred divergence times.
/// Why: Tree topology combined with temporal axis aids interpretation.
/// How: Horizontal layout with x-axis as calendar time, branches scaled by divergence times.
pub fn plot_time_tree(_graph: &GraphAncestral, _out_base: &Path) -> Result<(), Report> {
  todo!("Generate time-scaled tree visualization")
}
