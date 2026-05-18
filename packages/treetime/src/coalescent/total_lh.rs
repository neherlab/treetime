use crate::coalescent::edge_data::{collect_coalescent_edges, sum_coalescent_cost};
use crate::coalescent::events::collect_tree_events;
use crate::coalescent::integration::compute_integral_merger_rate;
use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::representation::payload::traits::TimetreeNode;
use eyre::Report;
use treetime_distribution::Distribution;
use treetime_graph::edge::{GraphEdge, TimeLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Computes the total coalescent log-likelihood of the tree under the given Tc.
///
/// For each non-root edge, evaluates the per-edge cost:
///   cost = I(t_merger) - I(t_node) - log(λ(t_merger)) * (m-1)/m
///
/// and returns LH = -Σ cost.
///
/// Multiplicity m is the number of children of the child node for internal nodes,
/// or 2.0 for leaves (binary merger assumption). This corrects a v0 erratum where
/// `total_LH()` uses fixed multiplicity=2 for all edges (see
/// `docs/port-v0-errata/coalescent-total-lh-fixed-multiplicity.md`).
///
/// This metric reports the full Kingman coalescent log-likelihood. The current
/// backward pass applies only internal-node contributions (leaf and root terms
/// are missing, see `docs/port-known-issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md`).
/// When those terms are added to inference, this metric will automatically align
/// with the optimized objective.
///
/// Accepts any `Distribution` for Tc (constant, skyline, or formula-based).
/// Piecewise functions are evaluated in TBP (time-before-present) coordinates.
pub fn compute_coalescent_total_lh<N, E, D>(graph: &Graph<N, E, D>, tc_dist: &Distribution) -> Result<f64, Report>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimeLength,
  D: Sync + Send,
{
  let (present_time, events_calendar) = collect_tree_events(graph)?;
  let events_tbp: Vec<_> = events_calendar
    .iter()
    .map(|(t, delta)| (t.to_tbp(present_time), *delta))
    .collect();

  let lineage_counts = compute_lineage_count_distribution(&events_tbp)?;
  let integral_merger_rate = compute_integral_merger_rate(tc_dist, &lineage_counts)?;
  let edges = collect_coalescent_edges(graph, present_time);

  sum_coalescent_cost(&edges, &integral_merger_rate, &lineage_counts, tc_dist)
}
