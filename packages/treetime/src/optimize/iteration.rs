use num_traits::pow::pow;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;

pub const DAMPING_FLOOR: f64 = 0.01;

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
