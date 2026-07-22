use eyre::Report;
use ndarray::Array1;
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;
use treetime_utils::make_error;

/// Per-branch coalescent merger rate κ(t) = (k(t)-1)/(2·Tc(t)).
///
/// This is the rate for one branch to merge with any other branch. The lineage
/// count is clamped as `max(0.5, k - 1)` (see [`compute_merger_rate_lineage_count`]).
pub(super) fn compute_merger_rate_per_lineage_scalar(k: f64, tc: f64) -> f64 {
  compute_merger_rate_per_lineage(compute_merger_rate_lineage_count(k), tc)
}

/// Total coalescent merger rate λ(t) = k(t)·(k(t)-1)/(2·Tc(t)).
///
/// This is the rate for any pair of branches to merge. The lineage count is
/// clamped as `max(0.5, k - 1)` (see [`compute_merger_rate_lineage_count`]).
pub(super) fn compute_merger_rate_total_scalar(k: f64, tc: f64) -> f64 {
  compute_merger_rate_total(compute_merger_rate_lineage_count(k), tc)
}

/// Computes H(t) = ∫ₜᴾ κ(t') dt' in decimal calendar years.
///
/// This integral represents the expected number of merger events experienced by a branch.
/// Uses the exact breakpoints from lineage_counts where the function has discontinuities.
///
/// The integral is zero at the most recent event P and increases into the past.
pub fn compute_integral_merger_rate(
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
) -> Result<PiecewiseLinearFn, Report> {
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.len() < 2 {
    return make_error!("lineage count must have at least 2 breakpoints");
  }

  let n = breakpoints.len();
  let mut integral_values = vec![0.0; n];

  for i in (0..n - 1).rev() {
    let t0 = breakpoints[i];
    let t1 = breakpoints[i + 1];
    let dt = t1 - t0;

    // k is constant between breakpoints, use value at midpoint
    let mid = f64::midpoint(t0, t1);
    let k = lineage_counts.eval(mid);
    let tc = tc_dist.eval(mid)?;

    if !k.is_finite() {
      return make_error!("Coalescent lineage count must be finite at calendar time {mid:.6e}, got {k:.6e}");
    }
    if !tc.is_finite() || tc <= 0.0 {
      return make_error!("Coalescent Tc must be finite and positive at calendar time {mid:.6e}, got {tc:.6e}");
    }
    let rate = compute_merger_rate_per_lineage_scalar(k, tc);
    integral_values[i] = integral_values[i + 1] + dt * rate;
  }

  Ok(PiecewiseLinearFn::new(
    Array1::from(breakpoints.to_vec()),
    Array1::from_vec(integral_values),
  ))
}

/// Computes the textbook per-lineage merger integral H₀(t) = ∫ₜᴾ (k(t')-1)/2 dt'.
///
/// This is the coalescent survival integrand with the time scale Tc factored out
/// (i.e. Tc = 1) and *without* the `max(0.5, k-1)` clamp applied in the likelihood
/// path. It underlies the analytic Tc optimum: the exponent of the coalescent
/// likelihood is `Σ_edges [H₀(parent) − H₀(child)] / Tc` (see
/// [`crate::coalescent::optimize_tc`]).
///
/// Like [`compute_integral_merger_rate`], the integral is zero at the most recent
/// event P and increases into the past, using the exact lineage-count breakpoints.
pub fn compute_bare_integral_merger_rate(lineage_counts: &PiecewiseConstantFn) -> Result<PiecewiseLinearFn, Report> {
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.len() < 2 {
    return make_error!("lineage count must have at least 2 breakpoints");
  }

  let n = breakpoints.len();
  let mut integral_values = vec![0.0; n];

  for i in (0..n - 1).rev() {
    let t0 = breakpoints[i];
    let t1 = breakpoints[i + 1];
    let dt = t1 - t0;

    // k is constant between breakpoints, use value at midpoint
    let mid = f64::midpoint(t0, t1);
    let k = lineage_counts.eval(mid);
    if !k.is_finite() {
      return make_error!("Coalescent lineage count must be finite at calendar time {mid:.6e}, got {k:.6e}");
    }
    let rate = 0.5 * (k - 1.0);
    integral_values[i] = integral_values[i + 1] + dt * rate;
  }

  Ok(PiecewiseLinearFn::new(
    Array1::from(breakpoints.to_vec()),
    Array1::from_vec(integral_values),
  ))
}

/// Clamped effective lineage count `max(0.5, k - 1)`, matching TreeTime v0.
///
/// At the present boundary, `k` can be zero before the first sample event; v0
/// documents that this clamped region is evaluated only when the tree changes.
/// NaN lineage counts propagate as they do through v0's `numpy.maximum` call.
fn compute_merger_rate_lineage_count(k: f64) -> f64 {
  let n_lineages = k - 1.0;
  if n_lineages.is_nan() {
    n_lineages
  } else {
    f64::max(0.5, n_lineages)
  }
}

fn compute_merger_rate_per_lineage(n_lineages: f64, tc: f64) -> f64 {
  0.5 * n_lineages / tc
}

fn compute_merger_rate_total(n_lineages: f64, tc: f64) -> f64 {
  0.5 * n_lineages * (n_lineages + 1.0) / tc
}
