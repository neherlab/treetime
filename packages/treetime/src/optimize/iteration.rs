use itertools::izip;
use num_traits::pow::pow;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;

/// Save current branch lengths for all edges in graph traversal order.
pub fn save_branch_lengths<N, E, D>(graph: &Graph<N, E, D>) -> Vec<f64>
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  graph
    .get_edges()
    .iter()
    .map(|edge_ref| edge_ref.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
    .collect::<Vec<_>>()
}

/// Restore previously saved branch lengths to all edges in graph traversal order.
///
/// Inverse of [`save_branch_lengths`]. The saved vector must have been produced
/// from the same graph (same edge count and order).
pub fn restore_branch_lengths<N, E, D>(graph: &Graph<N, E, D>, saved: &[f64])
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  debug_assert_eq!(
    graph.get_edges().len(),
    saved.len(),
    "restore_branch_lengths: edge count changed between save and restore"
  );
  for (edge_ref, &bl) in izip!(graph.get_edges(), saved) {
    edge_ref.write_arc().payload().write_arc().set_branch_length(Some(bl));
  }
}

/// Minimum fraction of the old branch length retained at any iteration.
///
/// Without a floor, exponential damping $d^{i+1}$ decays to effectively zero
/// at high iteration counts (e.g. $0.75^{20} \approx 0.003$). On datasets where
/// the sparse variable/fixed position boundary oscillates, fully undamped late
/// iterations amplify the discrete jump. The floor ensures at least 1% of the
/// old value is retained, bridging the discontinuity at all iteration counts.
pub const DAMPING_FLOOR: f64 = 0.01;

/// Blend optimized branch lengths with saved old values using exponential damping.
///
/// At iteration `i` (0-based), each branch length becomes:
///   bl = bl_optimized * (1 - damping_factor) + bl_old * damping_factor
///
/// where `damping_factor = max(damping^(i+1), DAMPING_FLOOR)`.
///
/// When `damping == 0.0`, damping_factor = 0 and the optimized value is kept unchanged.
/// Early iterations take conservative steps; later iterations approach the full update
/// but never go below the `DAMPING_FLOOR` weight on the old value.
pub fn apply_damping<N, E, D>(graph: &Graph<N, E, D>, old_branch_lengths: &[f64], damping: f64, iteration: usize)
where
  N: GraphNode,
  E: GraphEdge + HasBranchLength,
  D: Send + Sync,
{
  if damping == 0.0 {
    return;
  }
  debug_assert_eq!(
    graph.get_edges().len(),
    old_branch_lengths.len(),
    "apply_damping: edge count changed between save and apply"
  );
  let damping_factor = pow(damping, iteration + 1).max(DAMPING_FLOOR);
  let new_weight = 1.0 - damping_factor;
  for (edge_ref, &old_bl) in izip!(graph.get_edges(), old_branch_lengths) {
    let mut edge = edge_ref.write_arc().payload().write_arc();
    let optimized_bl = edge.branch_length().unwrap_or(0.0);
    let damped_bl = optimized_bl * new_weight + old_bl * damping_factor;
    edge.set_branch_length(Some(damped_bl));
  }
}
