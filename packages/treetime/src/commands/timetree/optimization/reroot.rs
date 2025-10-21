use crate::representation::partition_timetree::GraphTimetree;
use eyre::Report;

/// Find optimal root location for temporal signal maximization.
///
/// Optional: Can improve divergence time estimates when true root is unknown.
/// Why: Root-to-tip regression fit depends on root placement.
/// How: Test multiple roots, choose one maximizing R² in temporal regression.
pub fn reroot_tree(_graph: &GraphTimetree, _method: &str) -> Result<(), Report> {
  todo!("Rerooting: least-squares (minimize residuals), min_dev (minimize variation), oldest (root on oldest sample)")
}
