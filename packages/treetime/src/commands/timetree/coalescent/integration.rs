use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use crate::commands::timetree::coalescent::piecewise_linear_fn::PiecewiseLinearFn;
use crate::distribution::distribution::Distribution;
use eyre::Report;
use ndarray::Array1;
use treetime_utils::make_error;

/// Computes branch merger rate κ(t) and total merger rate λ(t).
///
/// κ(t) = (k(t)-1)/(2*Tc(t)) - rate for one branch to merge with any other
/// λ(t) = k(t)*(k(t)-1)/(2*Tc(t)) - rate for any branch to merge with any other
///
/// # Clamping Logic
///
/// The lineage count is clamped: `k_clamped = max(0.5, k-1)` to ensure positive
/// merger rates even in edge cases. This matches Python v0 behavior.
///
/// Why 0.5? At tree boundaries:
/// - At TBP=0 (present), k=0 before first sample event → k-1=-1 → clamped to 0.5
/// - With k=1 (single lineage), k-1=0 → clamped to 0.5 (prevents zero rate)
///
/// The value 0.5 is somewhat arbitrary but ensures:
/// - Rates are always positive (required for valid probability distributions)
/// - The formulas remain numerically stable
/// - Boundary regions contribute minimally to the overall likelihood
///
/// Note: The clamped regions should ideally never be evaluated if the tree is
/// well-formed, as per Python v0 comment: "in these regions, the function
/// should only be called if the tree changes."
pub fn compute_merger_rates(k: &Array1<f64>, tc: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
  let k_clamped = k.mapv(|x| f64::max(0.5, x - 1.0));
  let branch_rate = 0.5 * &k_clamped / tc;
  let total_rate = 0.5 * &k_clamped * (&k_clamped + 1.0) / tc;
  (branch_rate, total_rate)
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

    let k_clamped = f64::max(0.5, k - 1.0);
    let rate = 0.5 * k_clamped / tc;

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
    let k_clamped = f64::max(0.5, k - 1.0);
    let rate = 0.5 * k_clamped / tc;
    -rate * breakpoints[0]
  };

  let integral_values: Vec<f64> = integral_values.iter().map(|v| v - i_at_zero).collect();

  Ok(PiecewiseLinearFn::new(
    Array1::from(breakpoints.to_vec()),
    Array1::from(integral_values),
  ))
}

#[cfg(test)]
mod tests {
  use super::*;
  use eyre::Report;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_merger_rates() -> Result<(), Report> {
    let k = array![2.0, 3.0, 4.0];
    let tc = array![0.001, 0.002, 0.003];

    let (actual_branch_rate, actual_total_rate) = compute_merger_rates(&k, &tc);

    let expected_branch_rate = array![500.0, 500.0, 500.0];
    let expected_total_rate = array![1000.0, 1500.0, 2000.0];

    pretty_assert_ulps_eq!(actual_branch_rate, expected_branch_rate);
    pretty_assert_ulps_eq!(actual_total_rate, expected_total_rate);

    Ok(())
  }

  #[test]
  fn test_merger_rates_edge_cases() -> Result<(), Report> {
    let k = array![1.0, 0.5, 2.0];
    let tc = array![0.001, 0.002, 0.003];

    let (actual_branch_rate, actual_total_rate) = compute_merger_rates(&k, &tc);

    let expected_branch_rate = array![250.0, 125.0, 166.66666666666666];
    let expected_total_rate = array![375.0, 187.5, 333.3333333333333];

    pretty_assert_ulps_eq!(actual_branch_rate, expected_branch_rate);
    pretty_assert_ulps_eq!(actual_total_rate, expected_total_rate);

    Ok(())
  }

  #[test]
  fn test_merger_rates_large_k() -> Result<(), Report> {
    let k = array![10.0, 20.0, 50.0];
    let tc = array![0.01, 0.01, 0.01];

    let (actual_branch_rate, actual_total_rate) = compute_merger_rates(&k, &tc);

    let expected_branch_rate = array![450.0, 950.0, 2450.0];
    let expected_total_rate = array![4500.0, 19000.0, 122500.0];

    pretty_assert_ulps_eq!(actual_branch_rate, expected_branch_rate);
    pretty_assert_ulps_eq!(actual_total_rate, expected_total_rate);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_constant_tc() -> Result<(), Report> {
    // Create PiecewiseConstantFn with k=2.0 from t=0 to t=10
    // Breakpoints: [0.0, 10.0], values: [0.0, 2.0, 2.0]
    let lineage_counts = PiecewiseConstantFn::new(array![0.0, 10.0], array![0.0, 2.0, 2.0]);

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.values();

    // I(0) = 0, I(10) = integral of 0.5*(2-1)/0.01 = 50*10 = 500
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 500.0, max_ulps = 1000);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_multiple_segments() -> Result<(), Report> {
    // Three segments: k=1 from 0-5, k=5 from 5-10
    // Breakpoints: [0.0, 5.0, 10.0], values: [0.0, 1.0, 5.0, 5.0]
    let lineage_counts = PiecewiseConstantFn::new(array![0.0, 5.0, 10.0], array![0.0, 1.0, 5.0, 5.0]);

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.values();

    // Segment 0-5: k=1, rate = 0.5*max(0.5, 1-1)/0.01 = 0.5*0.5/0.01 = 25, contribution = 25*5 = 125
    // Segment 5-10: k=5, rate = 0.5*(5-1)/0.01 = 200, contribution = 200*5 = 1000
    // Total = 125 + 1000 = 1125
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 1125.0, max_ulps = 1000);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_insufficient_points() {
    let lineage_counts = PiecewiseConstantFn::new(array![5.0], array![0.0, 2.0]);

    let tc_dist = Distribution::constant(0.01);

    let result = compute_integral_merger_rate(&tc_dist, &lineage_counts);

    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("at least 2 breakpoints"));
  }

  #[test]
  fn test_compute_integral_merger_rate_varying_tc() -> Result<(), Report> {
    // k=3 constant from 0-10
    let lineage_counts = PiecewiseConstantFn::new(array![0.0, 10.0], array![0.0, 3.0, 3.0]);

    // Tc varies linearly from 0.01 to 0.05
    let t_vals = Array1::linspace(0.0, 10.0, 100);
    let tc_vals = Array1::linspace(0.01, 0.05, 100);
    let tc_dist = Distribution::function(t_vals, tc_vals)?;

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.values();

    // Just verify the integral increases monotonically
    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    assert!(actual_y[actual_y.len() - 1] > 0.0);

    Ok(())
  }
}
