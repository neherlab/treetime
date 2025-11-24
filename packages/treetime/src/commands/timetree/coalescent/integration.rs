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
  use eyre::Report;
  use ndarray::{Array1, array};
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
    let t_vals = Array1::linspace(0.0, 10.0, 100);
    let k_vals = Array1::ones(100) * 2.0;
    let lineage_counts = Distribution::function(t_vals, k_vals)?;

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.y();

    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 500.0, max_ulps = 1000);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_increasing_lineages() -> Result<(), Report> {
    let t_vals = Array1::linspace(0.0, 10.0, 100);
    let k_vals = Array1::linspace(1.0, 10.0, 100);
    let lineage_counts = Distribution::function(t_vals, k_vals)?;

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.y();

    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 2257.0018365472897);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_insufficient_points() {
    use ndarray::Array1;

    let t_vals = Array1::from_vec(vec![0.0]);
    let k_vals = Array1::from_vec(vec![2.0]);
    let lineage_counts = Distribution::function(t_vals, k_vals).unwrap();

    let tc_dist = Distribution::constant(0.01);

    let result = compute_integral_merger_rate(&tc_dist, &lineage_counts);

    assert!(result.is_err());
    assert!(result.unwrap_err().to_string().contains("at least 2 points"));
  }

  #[test]
  fn test_compute_integral_merger_rate_trapezoidal_accuracy() -> Result<(), Report> {
    let t_vals = Array1::linspace(0.0, 1.0, 1000);
    let k_vals = Array1::ones(1000) * 2.0;
    let lineage_counts = Distribution::function(t_vals, k_vals)?;

    let tc_dist = Distribution::constant(0.01);

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.y();

    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 50.0, max_ulps = 1000);

    Ok(())
  }

  #[test]
  fn test_compute_integral_merger_rate_varying_tc() -> Result<(), Report> {
    let t_vals = Array1::linspace(0.0, 10.0, 100);
    let k_vals = Array1::ones(100) * 3.0;
    let lineage_counts = Distribution::function(t_vals.clone(), k_vals)?;

    let tc_vals = Array1::linspace(0.01, 0.05, 100);
    let tc_dist = Distribution::function(t_vals, tc_vals)?;

    let actual = compute_integral_merger_rate(&tc_dist, &lineage_counts)?;

    let actual_y = actual.y();

    pretty_assert_ulps_eq!(actual_y[0], 0.0);
    pretty_assert_ulps_eq!(actual_y[actual_y.len() - 1], 402.3921222992276);

    Ok(())
  }
}
