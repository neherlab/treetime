use crate::representation::edge_timetree::EdgeTimetree;
use crate::representation::node_timetree::NodeTimetree;
use crate::representation::partition_timetree::{GraphTimetree, PartitionTreetimeMarginalOps};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;

/// Resolve multifurcations using temporal constraints and sequence data.
///
/// Optional: Improves tree topology when input tree has polytomies.
/// Why: Binary trees are more informative and required for some analyses.
/// How: Greedy optimization of likelihood gain from inserting internal nodes.
pub fn resolve_polytomies(
  _graph: &GraphTimetree,
  _partitions: &[Arc<RwLock<dyn PartitionTreetimeMarginalOps<NodeTimetree, EdgeTimetree>>>],
) -> Result<usize, Report> {
  todo!(
    "For each polytomy: test pairwise mergers, choose highest likelihood gain. Returns count of resolved polytomies."
  )
}

/// Reset internal state after tree topology changes.
///
/// Core: Required after polytomy resolution to maintain consistency.
/// Why: Topology changes invalidate cached likelihood values and node attributes.
/// How: Clear cached data, recalculate tree layout, reset branch flags.
pub fn prepare_tree_after_topology_change(_graph: &GraphTimetree) -> Result<(), Report> {
  todo!("Clear cached likelihoods, recalculate dist2root, reset bad_branch flags")
}
