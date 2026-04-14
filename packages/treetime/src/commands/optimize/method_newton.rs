use crate::commands::optimize::optimize_unified::{
  OptimizationContribution, OptimizationMetrics, evaluate_with_indels, grid_search_inner, newton_tolerance,
};
use num::clamp;

/// Relative tolerance for Newton inner-loop convergence.
pub(crate) const NEWTON_REL_TOL: f64 = 0.001;

/// Absolute tolerance floor for Newton inner-loop convergence (subs/site).
/// Prevents the purely relative tolerance from degenerating to zero when
/// the current branch length is zero or very small. The value 1e-8 is well
/// below the smallest meaningful branch length (1/L ~ 1e-4 for L=10000)
/// so it does not mask real convergence.
pub(crate) const NEWTON_ABS_TOL: f64 = 1e-8;

/// Maximum iterations for Newton inner loop.
pub(crate) const NEWTON_MAX_ITER: usize = 10;

/// Newton-Raphson inner loop in $t$-space.
///
/// Takes an initial `branch_length` and `metrics` (already evaluated at that point),
/// runs up to `NEWTON_MAX_ITER` Newton steps with the convergence criterion
/// $|t_{\text{new}} - t_{\text{old}}| < \text{tol}(t)$, and returns the optimized
/// branch length. Falls through to grid search if the Hessian becomes non-negative.
pub(crate) fn newton_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  if metrics.second_derivative < 0.0 {
    let mut bl = branch_length;
    let mut new_bl = (bl - clamp(metrics.derivative / metrics.second_derivative, -1.0, bl)).max(min_branch_length);

    for _ in 0..NEWTON_MAX_ITER {
      if (new_bl - bl).abs() <= newton_tolerance(bl) {
        break;
      }
      let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, new_bl);
      if new_metrics.second_derivative < 0.0 {
        bl = new_bl;
        new_bl = (bl - clamp(new_metrics.derivative / new_metrics.second_derivative, -1.0, bl)).max(min_branch_length);
      } else {
        break;
      }
    }
    new_bl
  } else {
    grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation)
  }
}

/// Newton-Raphson inner loop in $\sqrt{t}$ space.
///
/// Reparameterizes the optimization variable as $s = \sqrt{t}$ and applies the
/// chain rule to transform $t$-space derivatives into $s$-space:
///
///   $d\ell/ds = 2s \cdot d\ell/dt$
///   $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$
///
/// This reduces the indel Hessian singularity from $O(1/t^2)$ to $O(1/t)$,
/// allowing Newton to make progress toward the combined (substitution + indel)
/// optimum instead of stalling at the indel-only MLE.
///
/// The convergence criterion is applied in $s$-space:
/// $|s_{\text{new}} - s_{\text{old}}| < \text{tol}(s)$.
pub(crate) fn newton_sqrt_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  let mut s = branch_length.sqrt();
  let min_s = min_branch_length.sqrt();

  // Transform initial metrics to s-space
  let (ds, d2s) = chain_rule_sqrt(s, metrics.derivative, metrics.second_derivative);

  if d2s >= 0.0 {
    return grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation);
  }

  let mut new_s = (s - clamp(ds / d2s, -1.0, s)).max(min_s);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_s - s).abs() <= newton_tolerance(s) {
      break;
    }
    let t = new_s * new_s;
    let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, t);
    let (ds_new, d2s_new) = chain_rule_sqrt(new_s, new_metrics.derivative, new_metrics.second_derivative);
    if d2s_new < 0.0 {
      s = new_s;
      new_s = (s - clamp(ds_new / d2s_new, -1.0, s)).max(min_s);
    } else {
      break;
    }
  }

  new_s * new_s
}

/// Transform $t$-space derivatives to $\sqrt{t}$-space via the chain rule.
///
/// Given $s = \sqrt{t}$:
///   $d\ell/ds = 2s \cdot d\ell/dt$
///   $d^2\ell/ds^2 = 4s^2 \cdot d^2\ell/dt^2 + 2 \cdot d\ell/dt$
pub(crate) fn chain_rule_sqrt(s: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_ds = 2.0 * s * dl_dt;
  let d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
  (dl_ds, d2l_ds2)
}
