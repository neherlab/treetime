use crate::clock::find_best_root::cost_function::BranchPointCostFunction;
use crate::clock::find_best_root::find_best_root::FindRootResult;
use log::debug;
use treetime_graph::edge::GraphEdgeKey;

/// Below this curvature, relative to the scale of the scores, the score along the edge counts as linear
const FLAT_CURVATURE: f64 = 1e-12;

/// Root position on one edge that minimizes the score of an ordinary least-squares regression, in closed form.
///
/// Every tip has the same weight, and a root at fraction `x` from the parent shifts
/// the divergences of the tips below the edge by `(1 - x) l` and of all others by `x l`. The sums `sum d`, `sum t d`
/// are then linear and `sum d^2` quadratic in `x`, while `sum t`, `sum t^2` and the norm do not change, so both
/// objectives are quadratics, `score(x) = a x^2 + b x + c`, and the parabola through `x = 0, 1/2, 1` is exact:
/// `c = s(0)`, `a = 2 (s(0) - 2 s(1/2) + s(1))`, `b = s(1) - s(0) - a`. Its vertex `-b / 2a`, clamped to `[0, 1]`, is
/// the optimum. A required positive rate is checked on this optimum by the caller rather than imposed as a
/// constraint: a constrained optimum would sit where the rate crosses zero, so data without temporal signal would get
/// a rate near zero instead of the error that asks for a given rate.
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
