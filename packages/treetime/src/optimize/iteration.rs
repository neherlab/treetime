use num_traits::pow::pow;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;

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
///
/// Operates on the loop's branch-length map: `branch_lengths` holds the freshly optimized
/// lengths and is blended in place with the pre-optimization `old_branch_lengths`. Both maps are
/// keyed by the same edge set, so each optimized length is blended with its own old length. A
/// missing weight (`None`) resolves to `0.0` for the blend, matching the marginal-input derivation.
pub fn apply_damping(
  branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  old_branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  damping: f64,
  iteration: usize,
) {
  if damping == 0.0 {
    return;
  }
  let damping_factor = pow(damping, iteration + 1).max(DAMPING_FLOOR);
  let new_weight = 1.0 - damping_factor;
  for (key, bl) in branch_lengths.iter_mut() {
    let optimized_bl = bl.unwrap_or(0.0);
    let old_bl = old_branch_lengths[key].unwrap_or(0.0);
    *bl = Some(optimized_bl * new_weight + old_bl * damping_factor);
  }
}
