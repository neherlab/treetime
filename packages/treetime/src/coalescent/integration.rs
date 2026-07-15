use eyre::Report;
use ndarray::{Array1, Zip};
use treetime_distribution::Distribution;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_grid::piecewise_linear_fn::PiecewiseLinearFn;
use treetime_utils::make_error;

/// Per-lineage and total coalescent merger rates.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MergerRates<T> {
  /// Rate for one lineage to merge with any other lineage.
  pub per_lineage: T,
  /// Rate for any pair of lineages to merge.
  pub total: T,
}

/// Computes scalar per-lineage merger rate κ(t) and total merger rate λ(t).
///
/// κ(t) = (k(t)-1)/(2*Tc(t)) - rate for one branch to merge with any other
/// λ(t) = k(t)*(k(t)-1)/(2*Tc(t)) - rate for any branch to merge with any other
///
/// # Clamping Logic
///
/// The lineage count is clamped as `max(0.5, k - 1)`, matching TreeTime v0.
/// At the present boundary, `k` can be zero before the first sample event; v0
/// documents that this clamped region is evaluated only when the tree changes.
/// NaN lineage counts propagate as they do through v0's `numpy.maximum` call.
pub fn compute_merger_rates_scalar(k: f64, tc: f64) -> MergerRates<f64> {
  let n_lineages = compute_merger_rate_lineage_count(k);
  MergerRates {
    per_lineage: compute_merger_rate_per_lineage(n_lineages, tc),
    total: compute_merger_rate_total(n_lineages, tc),
  }
}

pub(super) fn compute_merger_rate_per_lineage_scalar(k: f64, tc: f64) -> f64 {
  compute_merger_rate_per_lineage(compute_merger_rate_lineage_count(k), tc)
}

pub(super) fn compute_merger_rate_total_scalar(k: f64, tc: f64) -> f64 {
  compute_merger_rate_total(compute_merger_rate_lineage_count(k), tc)
}

/// Computes array-valued per-lineage merger rate κ(t) and total merger rate λ(t).
///
/// Length-one inputs follow ndarray broadcasting semantics.
#[cfg_attr(not(test), allow(dead_code))]
pub fn compute_merger_rates(k: &Array1<f64>, tc: &Array1<f64>) -> MergerRates<Array1<f64>> {
  let output_len = if k.len() == 1 { tc.len() } else { k.len() };
  let mut per_lineage = Array1::zeros(output_len);
  let mut total = Array1::zeros(output_len);

  Zip::from(&mut per_lineage)
    .and(&mut total)
    .and_broadcast(k)
    .and_broadcast(tc)
    .for_each(|per_lineage, total, &k, &tc| {
      let rates = compute_merger_rates_scalar(k, tc);
      *per_lineage = rates.per_lineage;
      *total = rates.total;
    });

  MergerRates { per_lineage, total }
}

/// Computes I(t) = ∫₀ᵗ κ(t') dt' via piecewise integration.
///
/// This integral represents the expected number of merger events experienced by a branch.
/// Uses the exact breakpoints from lineage_counts where the function has discontinuities.
///
/// The integral is normalized so that I(0) = 0, meaning the integral starts at the present
/// (TBP = 0) and increases going into the past.
pub fn compute_integral_merger_rate(
  tc_dist: &Distribution,
  lineage_counts: &PiecewiseConstantFn,
) -> Result<PiecewiseLinearFn, Report> {
  let breakpoints = lineage_counts.breakpoints();
  if breakpoints.len() < 2 {
    return make_error!("lineage count must have at least 2 breakpoints");
  }

  let n = breakpoints.len();
  let mut integral_values = Vec::with_capacity(n);
  integral_values.push(0.0); // I(breakpoints[0]) = 0

  for i in 1..n {
    let t0 = breakpoints[i - 1];
    let t1 = breakpoints[i];
    let dt = t1 - t0;

    // k is constant between breakpoints, use value at midpoint
    let mid = f64::midpoint(t0, t1);
    let k = lineage_counts.eval(mid);
    let tc = tc_dist.eval(mid)?;

    let rate = compute_merger_rate_per_lineage_scalar(k, tc);

    integral_values.push(integral_values[i - 1] + dt * rate);
  }

  // First breakpoint is at earliest sample (TBP ~ 0)
  // Normalize so I(0) = 0 by subtracting I(breakpoints[0])
  // For typical trees, breakpoints[0] is very close to 0
  // If breakpoints[0] > 0, we need to find I(0) by extrapolation
  let i_at_zero = if breakpoints[0] <= 0.0 {
    integral_values[0] // which is 0
  } else {
    // Extrapolate backwards: I(0) = I(bp[0]) - rate_0 * bp[0]
    let k = lineage_counts.eval(breakpoints[0] / 2.0);
    let tc = tc_dist.eval(breakpoints[0] / 2.0)?;
    let rate = compute_merger_rate_per_lineage_scalar(k, tc);
    -rate * breakpoints[0]
  };

  let integral_values: Vec<f64> = integral_values.iter().map(|v| v - i_at_zero).collect();

  Ok(PiecewiseLinearFn::new(
    Array1::from(breakpoints.to_vec()),
    Array1::from(integral_values),
  ))
}

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
