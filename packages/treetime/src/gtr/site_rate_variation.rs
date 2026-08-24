#[cfg(test)]
mod __tests__;

use eyre::Report;
use ndarray::Array1;
use statrs::distribution::{ContinuousCDF, Gamma};
use treetime_utils::{make_error, make_report};

/// Compute K discrete rate categories approximating a Gamma(alpha, alpha) distribution.
///
/// Returns K rate multipliers with mean 1.0 that approximate continuous gamma-distributed
/// rate variation across alignment sites. Each category represents the conditional mean
/// rate within one of K equal-probability quantile intervals.
///
/// # Background (Yang 1994)
///
/// Among-site rate variation (ASRV) models the observation that different alignment
/// positions evolve at different rates due to varying selective constraints.
/// The discrete gamma approximation divides the continuous Gamma(alpha, alpha)
/// distribution into K equal-probability intervals and represents each by its
/// conditional mean. This is the "+Γ" suffix in model notation (e.g., "GTR+Γ4").
///
/// The shape parameter alpha controls the degree of rate variation:
/// - alpha < 1: strong rate heterogeneity (many slow sites, few fast sites)
/// - alpha = 1: exponential distribution
/// - alpha > 1: moderate variation, bell-shaped
/// - alpha → ∞: uniform rates (no variation)
///
/// Values between 0.1 and 2.0 are typical for biological sequences.
///
/// # Algorithm
///
/// For Gamma(alpha, alpha) with mean = 1:
/// 1. Divide into K equal-probability intervals [q_{k-1}, q_k]
/// 2. Rate for category k = K * integral_{q_{k-1}}^{q_k} x * f(x) dx
///    = K * (F_{alpha+1,alpha}(q_k) - F_{alpha+1,alpha}(q_{k-1}))
///    where F_{alpha+1,alpha} is the CDF of Gamma(alpha+1, alpha).
///
/// # Reference
///
/// Yang Z (1994). "Maximum likelihood phylogenetic estimation from DNA sequences
/// with variable rates over sites: approximate methods." J Mol Evol 39:306-314.
/// DOI: 10.1007/BF00160154
pub fn discrete_gamma_rates(alpha: f64, n_categories: usize) -> Result<Array1<f64>, Report> {
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

  // Gamma(alpha, alpha) has mean = alpha/alpha = 1.0
  let gamma =
    Gamma::new(alpha, alpha).map_err(|e| make_report!("Failed to create Gamma({alpha}, {alpha}) distribution: {e}"))?;

  // Gamma(alpha+1, alpha) for computing conditional means within quantile intervals.
  // Derivation: integral of x * f_{a,b}(x) dx = (a/b) * F_{a+1,b}(x)
  // For a=alpha, b=alpha: integral = 1 * F_{alpha+1,alpha}(x)
  let gamma_next = Gamma::new(alpha + 1.0, alpha)
    .map_err(|e| make_report!("Failed to create Gamma({}, {alpha}) distribution: {e}", alpha + 1.0))?;

  let k = n_categories as f64;
  let mut rates = Array1::zeros(n_categories);

  for i in 0..n_categories {
    // Quantile boundaries for equal-probability intervals
    let q_lower = if i == 0 { 0.0 } else { gamma.inverse_cdf(i as f64 / k) };
    let q_upper = if i == n_categories - 1 {
      f64::INFINITY
    } else {
      gamma.inverse_cdf((i + 1) as f64 / k)
    };

    // Conditional mean in [q_lower, q_upper]:
    // r_k = K * (F_{alpha+1,alpha}(q_upper) - F_{alpha+1,alpha}(q_lower))
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
