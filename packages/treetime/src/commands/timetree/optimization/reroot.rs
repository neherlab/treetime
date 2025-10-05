use crate::commands::timetree::data::date_constraints::DateConstraintSet;
use crate::representation::graph_ancestral::GraphAncestral;
use eyre::Report;

/// Reroot tree to optimize temporal signal using root-to-tip regression.
///
/// Core: Required for molecular clock analysis when root position affects temporal signal quality.
/// Why: Finding optimal root maximizes correlation between sampling dates and genetic divergence.
/// How: Tests all possible root positions, evaluates regression fit for each.
pub fn reroot_tree(_graph: &GraphAncestral, _constraints: &DateConstraintSet, _method: &str) -> Result<(), Report> {
  todo!("Rerooting: least-squares (minimize residuals), min_dev (minimize variation), oldest (root on oldest sample)")
}
