use crate::commands::timetree::coalescent::events::collect_tree_events;
use crate::commands::timetree::coalescent::integration::compute_integral_merger_rate;
use crate::commands::timetree::coalescent::lineage_dynamics::compute_lineage_count_distribution;
use crate::commands::timetree::coalescent::time_coordinate::{CalendarTime, Tbp};
use crate::commands::timetree::timetree_traits::TimetreeNode;
use eyre::Report;
use log::warn;
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
/// or 2.0 for leaves (binary merger assumption).
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

  // Collect per-edge data first (graph traversal closure cannot return Result)
  let mut edges: Vec<EdgeBranchData> = Vec::new();

  graph.iter_breadth_first_forward(|node| {
    if node.parent_keys.is_empty() {
      return;
    }

    let Some(time_dist) = node.payload.time_distribution() else {
      return;
    };

    let Some(t_calendar) = time_dist.likely_time() else {
      return;
    };

    let t_node = CalendarTime::new(t_calendar).to_tbp(present_time);

    let (parent_node_key, parent_edge_key) = node.parent_keys[0];

    let branch_length = graph
      .get_edge(parent_edge_key)
      .and_then(|e| e.read_arc().payload().read_arc().time_length())
      .or_else(|| {
        graph.get_node(parent_node_key).and_then(|parent| {
          parent
            .read_arc()
            .payload()
            .read_arc()
            .time_distribution()
            .as_ref()
            .and_then(|d| d.likely_time())
            .map(|parent_t| {
              let parent_tbp = CalendarTime::new(parent_t).to_tbp(present_time);
              parent_tbp - t_node
            })
        })
      });

    let Some(branch_length) = branch_length else {
      warn!(
        "Coalescent total LH: skipping node (key={:?}) with undetermined branch length",
        node.key
      );
      return;
    };

    let multiplicity = if node.child_edges.is_empty() {
      2.0
    } else {
      node.child_edges.len() as f64
    };

    edges.push(EdgeBranchData {
      t_node,
      branch_length,
      multiplicity,
    });
  });

  // Compute total LH from collected edge data
  let mut total_lh = 0.0;

  for edge in &edges {
    let t_merger = edge.t_node + edge.branch_length.max(0.0);

    let i_merger = integral_merger_rate.eval(t_merger.value());
    let i_node = integral_merger_rate.eval(edge.t_node.value());

    let k = lineage_counts.eval(t_merger.value());
    let k_clamped = f64::max(0.5, k - 1.0);

    let tc_val = tc_dist.eval(t_merger.value())?;
    let lambda = 0.5 * k_clamped * (k_clamped + 1.0) / tc_val;

    let cost = (i_merger - i_node) - lambda.ln() * (edge.multiplicity - 1.0) / edge.multiplicity;
    total_lh -= cost;
  }

  Ok(total_lh)
}

struct EdgeBranchData {
  t_node: Tbp,
  branch_length: f64,
  multiplicity: f64,
}
