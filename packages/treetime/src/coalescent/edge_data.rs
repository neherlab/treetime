use crate::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use crate::coalescent::piecewise_linear_fn::PiecewiseLinearFn;
use crate::coalescent::time_coordinate::{CalendarTime, Tbp};
use crate::payload::traits::TimetreeNode;
use eyre::Report;
use log::warn;
use treetime_distribution::Distribution;
use treetime_graph::edge::{GraphEdge, TimeLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::make_error;

/// Per-edge data needed for coalescent likelihood computation.
///
/// Collected once from the graph, then reused for cost evaluation
/// (many times per Tc optimization run via Brent's method).
pub struct CoalescentEdgeData {
  /// Time before present of the child node.
  pub t_node: Tbp,
  /// Branch length to parent (time units, TBP direction).
  pub branch_length: f64,
  /// Number of children of the child node (2.0 for leaves, actual count for internal).
  pub multiplicity: f64,
}

/// Collects per-edge coalescent data from all non-root nodes.
///
/// For each non-root node, extracts the TBP time, branch length (from edge
/// time_length or parent-child time difference), and multiplicity (actual
/// child count for internal nodes, 2.0 for leaves under the binary merger
/// assumption).
pub fn collect_coalescent_edges<N, E, D>(graph: &Graph<N, E, D>, present_time: CalendarTime) -> Vec<CoalescentEdgeData>
where
  N: GraphNode + TimetreeNode,
  E: GraphEdge + TimeLength,
  D: Sync + Send,
{
  let mut edges = Vec::new();

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
        "Coalescent edge data: skipping node (key={:?}) with undetermined branch length",
        node.key
      );
      return;
    };

    // Multiplicity = number of children at the PARENT (merger) node.
    // The merger rate credit log(λ) * (m-1)/m is distributed across the m edges
    // entering the parent. Each edge gets the same share based on the parent's
    // child count, regardless of whether the child is a leaf or internal node.
    let multiplicity = graph
      .get_node(parent_node_key)
      .map_or(2.0, |parent| parent.read_arc().outbound().len() as f64);

    edges.push(CoalescentEdgeData {
      t_node,
      branch_length,
      multiplicity,
    });
  });

  edges
}

/// Sums per-edge coalescent cost to produce the total log-likelihood.
///
/// For each edge, the cost is:
///   cost = I(t_merger) - I(t_node) - log(λ(t_merger)) * (m-1)/m
///
/// Returns LH = -Σ cost.
///
/// `tc_dist` is evaluated at each merger time to support time-varying Tc
/// (skyline distributions). For constant Tc, `Distribution::constant(tc)`
/// returns the amplitude directly.
pub fn sum_coalescent_cost(
  edges: &[CoalescentEdgeData],
  integral_merger_rate: &PiecewiseLinearFn,
  lineage_counts: &PiecewiseConstantFn,
  tc_dist: &Distribution,
) -> Result<f64, Report> {
  let mut total_lh = 0.0;

  for edge in edges {
    let t_merger = edge.t_node + edge.branch_length.max(0.0);

    let i_merger = integral_merger_rate.eval(t_merger.value());
    let i_node = integral_merger_rate.eval(edge.t_node.value());

    // Use eval_left to get the pre-event lineage count at merger times.
    // Merger times fall exactly on PiecewiseConstantFn breakpoints (t_merger = parent_tbp).
    // The merger rate λ depends on lineages BEFORE the merger event.
    let k = lineage_counts.eval_left(t_merger.value());
    let k_clamped = f64::max(0.5, k - 1.0);

    let tc_val = tc_dist.eval(t_merger.value())?;
    if tc_val <= 0.0 {
      return make_error!(
        "Coalescent Tc must be positive, got {tc_val:.6e} at t={:.6e}",
        t_merger.value()
      );
    }
    let lambda = 0.5 * k_clamped * (k_clamped + 1.0) / tc_val;

    let cost = (i_merger - i_node) - lambda.ln() * (edge.multiplicity - 1.0) / edge.multiplicity;
    total_lh -= cost;
  }

  Ok(total_lh)
}
