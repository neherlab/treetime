use crate::distribution::distribution::Distribution;
use eyre::Report;
use ndarray::Array1;
use treetime_utils::make_error;

/// Computes branch merger rate κ(t) and total merger rate λ(t).
///
/// κ(t) = (k(t)-1)/(2*Tc(t)) - rate for one branch to merge with any other
/// λ(t) = k(t)*(k(t)-1)/(2*Tc(t)) - rate for any branch to merge with any other
///
/// k is clamped to ensure positive rates even in edge cases.
pub fn compute_merger_rates(k: &Array1<f64>, tc: &Array1<f64>) -> (Array1<f64>, Array1<f64>) {
  let k_clamped = k.mapv(|x| f64::max(0.5, x - 1.0));
  let branch_rate = 0.5 * &k_clamped / tc;
  let total_rate = 0.5 * &k_clamped * (&k_clamped + 1.0) / tc;
  (branch_rate, total_rate)
}

/// Computes I(t) = ∫₀ᵗ κ(t') dt' via trapezoidal integration.
///
/// This integral represents the expected number of merger events experienced by a branch.
/// Uses the exact time points from lineage_counts where the function has discontinuities.
pub fn compute_integral_merger_rate(
  tc_dist: &Distribution,
  lineage_counts: &Distribution,
) -> Result<Distribution, Report> {
  let tvals = lineage_counts.t();
  if tvals.len() < 2 {
    return make_error!("lineage count distribution must have at least 2 points");
  }

  let k_vals = lineage_counts.eval_many(&tvals)?;
  let tc_vals = tc_dist.eval_many(&tvals)?;
  let k_clamped = k_vals.mapv(|x| f64::max(0.5, x - 1.0));
  let branch_rates = 0.5 * &k_clamped / &tc_vals;

  let n_points = tvals.len();
  let mut integral_values = Array1::zeros(n_points);
  for i in 1..n_points {
    let dt = tvals[i] - tvals[i - 1];
    let avg_rate = 0.5 * (branch_rates[i - 1] + branch_rates[i]);
    integral_values[i] = integral_values[i - 1] + dt * avg_rate;
  }

  Distribution::function(tvals.to_owned(), integral_values)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_merger_rates() -> Result<(), Report> {
    let k = array![2.0, 3.0, 4.0];
    let tc = array![0.001, 0.002, 0.003];

    let (branch_rate, total_rate) = compute_merger_rates(&k, &tc);

    assert_abs_diff_eq!(branch_rate[0], 0.5 * 1.0 / 0.001, epsilon = 1e-10);
    assert_abs_diff_eq!(branch_rate[1], 0.5 * 2.0 / 0.002, epsilon = 1e-10);
    assert_abs_diff_eq!(branch_rate[2], 0.5 * 3.0 / 0.003, epsilon = 1e-10);

    assert_abs_diff_eq!(total_rate[0], 0.5 * 1.0 * 2.0 / 0.001, epsilon = 1e-10);
    assert_abs_diff_eq!(total_rate[1], 0.5 * 2.0 * 3.0 / 0.002, epsilon = 1e-10);
    assert_abs_diff_eq!(total_rate[2], 0.5 * 3.0 * 4.0 / 0.003, epsilon = 1e-10);

    Ok(())
  }
}
