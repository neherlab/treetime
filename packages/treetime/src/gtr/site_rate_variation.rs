#[cfg(test)]
mod __tests__;

use eyre::Report;
use ndarray::Array1;
use statrs::distribution::{ContinuousCDF, Gamma};
use treetime_utils::{make_error, make_report};

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn discrete_gamma_rates(alpha: f64, n_categories: usize) -> Result<Array1<f64>, Report> {
  if alpha < 0.15 {
    return make_error!(
      "Gamma shape parameter alpha must be >= 0.15 (statrs Gamma CDF is numerically \
       unstable for smaller values), got {alpha}"
    );
  }
  if n_categories == 0 {
    return make_error!("Number of rate categories must be at least 1, got {n_categories}");
  }
  if n_categories == 1 {
    return Ok(Array1::ones(1));
  }

  let gamma =
    Gamma::new(alpha, alpha).map_err(|e| make_report!("Failed to create Gamma({alpha}, {alpha}) distribution: {e}"))?;

  let gamma_next = Gamma::new(alpha + 1.0, alpha)
    .map_err(|e| make_report!("Failed to create Gamma({}, {alpha}) distribution: {e}", alpha + 1.0))?;

  let k = n_categories as f64;
  let mut rates = Array1::zeros(n_categories);

  for i in 0..n_categories {
    let q_lower = if i == 0 { 0.0 } else { gamma.inverse_cdf(i as f64 / k) };
    let q_upper = if i == n_categories - 1 {
      f64::INFINITY
    } else {
      gamma.inverse_cdf((i + 1) as f64 / k)
    };

    let cdf_lower = if q_lower == 0.0 { 0.0 } else { gamma_next.cdf(q_lower) };
    let cdf_upper = if q_upper == f64::INFINITY {
      1.0
    } else {
      gamma_next.cdf(q_upper)
    };

    rates[i] = k * (cdf_upper - cdf_lower);
  }

  Ok(rates)
}
