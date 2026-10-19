use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_root::FindRootResult;
use log::debug;
use treetime_graph::edge::GraphEdgeKey;

const FLAT_CURVATURE: f64 = 1e-12;

pub(crate) fn optimize_exact(edge: GraphEdgeKey, cost_fn: &BranchPointCostFunction) -> FindRootResult {
  let score_at = |x: f64| cost_fn.score_clock_set(&cost_fn.evaluate_clock_set(x));
  let (s0, s_half, s1) = (score_at(0.0), score_at(0.5), score_at(1.0));
  let a = 2.0 * (s0 - 2.0 * s_half + s1);
  let b = s1 - s0 - a;
  let scale = s0.abs().max(s1.abs()).max(f64::MIN_POSITIVE);
  let split = if a > FLAT_CURVATURE * scale {
    (-b / (2.0 * a)).clamp(0.0, 1.0)
  } else if s1 < s0 {
    1.0
  } else {
    0.0
  };
  let clock_set = cost_fn.evaluate_clock_set(split);
  let chisq = cost_fn.score_clock_set(&clock_set);
  debug!("Exact root position on edge {edge}: split = {split:.6}, chi-squared = {chisq:.6e}");
  FindRootResult {
    edge: Some(edge),
    split,
    clock_set,
    chisq,
  }
}
