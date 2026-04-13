use crate::commands::optimize::optimize_unified::{
  OptimizationContribution, OptimizationMetrics, evaluate_with_indels, grid_search_inner, newton_tolerance_log,
  newton_tolerance_sqrt, newton_tolerance_t,
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
///
/// **Convergence tolerance:** Step-size in $t$-space. Relative 0.1% of current $t$
/// with absolute floor 1e-8 subs/site: $\text{tol}(t) = \max(0.001 \cdot t, 10^{-8})$.
/// The tolerance is the same in both parameterized space and $t$-space (identity map).
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
      if (new_bl - bl).abs() <= newton_tolerance_t(bl) {
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
/// **Convergence tolerance:** Step-size in $s$-space: $\text{tol}(s) = \max(0.0005 \cdot s, 10^{-8})$.
/// Maps to tighter $t$-tolerance near zero because $ds/dt = 1/(2\sqrt{t})$ amplifies
/// precision: a step $\Delta s$ in $s$-space corresponds to $\Delta t \approx 2s \cdot \Delta s$
/// in $t$-space. At $s = 0.1$ ($t = 0.01$), tolerance 5e-5 in $s$ yields ~1e-5 in $t$.
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

  let mut new_s = (s - clamp(ds / d2s, sqrt_step_lower_bound(s), s)).max(min_s);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_s - s).abs() <= newton_tolerance_sqrt(s) {
      break;
    }
    let t = new_s * new_s;
    let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, t);
    let (ds_new, d2s_new) = chain_rule_sqrt(new_s, new_metrics.derivative, new_metrics.second_derivative);
    if d2s_new < 0.0 {
      s = new_s;
      new_s = (s - clamp(ds_new / d2s_new, sqrt_step_lower_bound(s), s)).max(min_s);
    } else {
      break;
    }
  }

  new_s * new_s
}

/// Lower bound on the Newton step $\delta_s$ in $\sqrt{t}$-space.
///
/// Limits $t$-space increase to 1.0 subs/site per iteration.
/// Derived from $(s - \delta_s)^2 - s^2 \leq 1$, giving
/// $\delta_s \geq s - \sqrt{s^2 + 1}$.
/// Always negative (permits $s$ to increase).
pub(crate) fn sqrt_step_lower_bound(s: f64) -> f64 {
  s - (s * s + 1.0).sqrt()
}

/// Lower bound on the Newton step $\delta_u$ in $\ln(t)$-space.
///
/// Limits $t$-space increase to 1.0 subs/site per iteration.
/// Derived from $e^{u - \delta_u} - e^u \leq 1$, giving
/// $\delta_u \geq -\ln(1 + 1/t)$. Uses `ln_1p` for accuracy at small $1/t$.
/// Always negative for $t > 0$ (permits $u$ to increase).
pub(crate) fn log_step_lower_bound(t: f64) -> f64 {
  -(1.0 / t).ln_1p()
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

/// Transform $t$-space derivatives to $\ln(t)$-space via the chain rule.
///
/// Given $u = \ln(t)$, the chain rule with $dt/du = t$ gives:
///   $d\ell/du = t \cdot d\ell/dt$
///   $d^2\ell/du^2 = t^2 \cdot d^2\ell/dt^2 + t \cdot d\ell/dt$
///
/// The Poisson indel Hessian in $u$-space becomes $-\mu t$ (bounded),
/// eliminating the $O(1/t^2)$ singularity that dominates in $t$-space.
pub(crate) fn chain_rule_log(t: f64, dl_dt: f64, d2l_dt2: f64) -> (f64, f64) {
  let dl_du = t * dl_dt;
  let d2l_du2 = t * t * d2l_dt2 + t * dl_dt;
  (dl_du, d2l_du2)
}

/// Newton-Raphson inner loop in $\ln(t)$ space.
///
/// Reparameterizes the optimization variable as $u = \ln(t)$ and applies the
/// chain rule to transform $t$-space derivatives into $u$-space:
///
///   $d\ell/du = t \cdot d\ell/dt$
///   $d^2\ell/du^2 = t^2 \cdot d^2\ell/dt^2 + t \cdot d\ell/dt$
///
/// This eliminates the indel Hessian singularity entirely: the Poisson indel
/// Hessian in $u$-space is $-\mu t$ (bounded for all $t > 0$), compared to
/// $-k/t^2$ in $t$-space and $-k/t$ in $\sqrt{t}$-space.
///
/// Step clamping in $u$-space:
/// - Upper bound $u - u_{\min}$: prevents $u_{\text{new}} < u_{\min}$
/// - Lower bound $-\ln(1 + 1/t)$: limits $t$-space increase to 1.0 subs/site
///
/// **Convergence tolerance:** Step-size in $u$-space: $\text{tol}(u) = \ln(1 + 0.001) \approx 0.001$.
/// This is a natural relative tolerance in $t$-space because $dt/t \approx du$ for small steps.
/// A step $\Delta u$ in $u$-space corresponds to $\Delta t \approx t \cdot \Delta u$ in $t$-space,
/// so the same tolerance achieves ~0.1% relative precision regardless of $t$.
pub(crate) fn newton_log_inner(
  branch_length: f64,
  metrics: &OptimizationMetrics,
  contributions: &[OptimizationContribution],
  indel_count: usize,
  indel_rate: f64,
  min_branch_length: f64,
  one_mutation: f64,
) -> f64 {
  debug_assert!(
    branch_length > 0.0,
    "newton_log_inner requires branch_length > 0, got {branch_length}"
  );

  let mut u = branch_length.ln();
  let u_min = min_branch_length.max(1e-12).ln();

  // Transform initial metrics to u-space
  let (du, d2u) = chain_rule_log(branch_length, metrics.derivative, metrics.second_derivative);

  if d2u >= 0.0 {
    return grid_search_inner(branch_length, contributions, indel_count, indel_rate, one_mutation);
  }

  // Step clamping bounds in u-space:
  // Lower: log_step_lower_bound(t) limits t-space increase to 1.0 subs/site
  // Upper: u - u_min prevents u from going below u_min
  let step_lower = log_step_lower_bound(branch_length);
  let step_upper = u - u_min;
  let mut new_u = (u - clamp(du / d2u, step_lower, step_upper)).max(u_min);

  for _ in 0..NEWTON_MAX_ITER {
    if (new_u - u).abs() <= newton_tolerance_log() {
      break;
    }
    let t = new_u.exp();
    let new_metrics = evaluate_with_indels(contributions, indel_count, indel_rate, t);
    let (du_new, d2u_new) = chain_rule_log(t, new_metrics.derivative, new_metrics.second_derivative);
    if d2u_new < 0.0 {
      u = new_u;
      let step_lower = log_step_lower_bound(t);
      let step_upper = u - u_min;
      new_u = (u - clamp(du_new / d2u_new, step_lower, step_upper)).max(u_min);
    } else {
      break;
    }
  }

  new_u.exp()
}
