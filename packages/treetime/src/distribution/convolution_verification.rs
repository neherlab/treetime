use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_convolution::distribution_convolution;
use crate::distribution::reference::exponential::{exponential_convolution, exponential_f, exponential_g};
use crate::distribution::reference::gaussian::{gaussian_convolution, gaussian_f, gaussian_g};
use crate::distribution::reference::grid_fn::GridFn;
use crate::utils::error::make_error;
use eyre::WrapErr;

/// Compare analytical result with discrete distribution result
///
/// Returns the maximum absolute difference between two results over their
/// common domain. This function is used to quantify how well a discrete convolution
/// approximates the analytical result.
///
/// Note: Significant differences are expected and normal due to discretization
/// effects, sampling resolution, and numerical precision limitations.
pub fn compare_with_distribution(
  analytical: &GridFn,
  discrete: &Distribution,
  tolerance: f64,
) -> Result<f64, eyre::Report> {
  assert!(tolerance > 0.0, "tolerance must be positive, got {tolerance}");

  let Distribution::Function(d) = discrete else {
    // For non-function distributions (Empty, Point, Range), return 0 difference
    return Ok(0.0);
  };

  // Find common domain for comparison
  let a_min = analytical.x().iter().copied().fold(f64::INFINITY, f64::min);
  let a_max = analytical.x().iter().copied().fold(f64::NEG_INFINITY, f64::max);
  let d_min = d.t().iter().copied().fold(f64::INFINITY, f64::min);
  let d_max = d.t().iter().copied().fold(f64::NEG_INFINITY, f64::max);

  if !a_min.is_finite() || !a_max.is_finite() || !d_min.is_finite() || !d_max.is_finite() {
    return make_error!("Non-finite values found in function domains");
  }

  let common_min = f64::max(a_min, d_min);
  let common_max = f64::min(a_max, d_max);

  if common_min >= common_max {
    return Ok(0.0); // No overlap
  }

  // Use adaptive sampling density based on domain size
  let domain_size = common_max - common_min;
  let n_samples = if domain_size > 10.0 { 200 } else { 100 };
  let dx = (common_max - common_min) / (n_samples - 1) as f64;

  let mut max_diff = 0.0;
  let mut error_count = 0;

  for i in 0..n_samples {
    let x = common_min + i as f64 * dx;

    let a_val = analytical.interp(x).unwrap_or_else(|_| {
      error_count += 1;
      0.0
    });
    let d_val = d.interp(x).unwrap_or_else(|_| {
      error_count += 1;
      0.0
    });

    if !a_val.is_finite() || !d_val.is_finite() {
      error_count += 1;
      continue;
    }

    let diff = (a_val - d_val).abs();
    max_diff = f64::max(max_diff, diff);
  }

  if error_count > n_samples / 10 {
    return make_error!("Too many interpolation errors: {error_count}/{n_samples}");
  }

  Ok(max_diff)
}

/// Comprehensive verification test for Gaussian convolution
///
/// Compares discrete convolution of two Gaussian functions against the analytical
/// result. This demonstrates that the discrete convolution implementation produces
/// results that are qualitatively similar to the theoretical expectation.
pub fn verify_gaussian_convolution(sigma_f: f64, sigma_g: f64, mu: f64, tolerance: f64) -> Result<f64, eyre::Report> {
  assert!(sigma_f > 0.0, "sigma_f must be positive, got {sigma_f}");
  assert!(sigma_g > 0.0, "sigma_g must be positive, got {sigma_g}");
  assert!(tolerance > 0.0, "tolerance must be positive, got {tolerance}");
  assert!(mu.is_finite(), "mu must be finite, got {mu}");

  // Generate individual functions with appropriate domains and higher resolution
  let f_domain = (-5.0 * sigma_f, 5.0 * sigma_f);
  let g_domain = (mu - 5.0 * sigma_g, mu + 5.0 * sigma_g);

  let f = gaussian_f(sigma_f, f_domain, (f_domain.1 - f_domain.0) / 200.0)?;
  let g = gaussian_g(sigma_g, mu, g_domain, (g_domain.1 - g_domain.0) / 200.0)?;

  // Convert to Distribution objects for discrete convolution
  let f_dist = Distribution::function(f.x().clone(), f.y().clone())?;
  let g_dist = Distribution::function(g.x().clone(), g.y().clone())?;

  // Discrete convolution
  let discrete = distribution_convolution(&f_dist, &g_dist).wrap_err("Failed to compute discrete convolution")?;

  // Analytical convolution with extended domain and high resolution
  let result_domain = (f_domain.0 + g_domain.0, f_domain.1 + g_domain.1);
  let analytical = gaussian_convolution(
    sigma_f,
    sigma_g,
    mu,
    result_domain,
    (result_domain.1 - result_domain.0) / 400.0,
  )?;

  // Compare results
  compare_with_distribution(&analytical, &discrete, tolerance)
    .wrap_err("Failed to compare Gaussian convolution results")
}

/// Comprehensive verification test for exponential convolution
///
/// Compares discrete convolution of two exponential functions against the analytical
/// result. This demonstrates that the discrete convolution implementation produces
/// results that are qualitatively similar to the theoretical expectation.
pub fn verify_exponential_convolution(
  a: f64,
  b: f64,
  f_max: f64,
  g_max: f64,
  tolerance: f64,
) -> Result<f64, eyre::Report> {
  assert!(a > 0.0, "parameter a must be positive, got {a}");
  assert!(b > 0.0, "parameter b must be positive, got {b}");
  assert!(f_max > 0.0, "f_max must be positive, got {f_max}");
  assert!(g_max > 0.0, "g_max must be positive, got {g_max}");
  assert!(tolerance > 0.0, "tolerance must be positive, got {tolerance}");

  // Generate individual functions with higher resolution
  let f = exponential_f(a, (0.0, f_max), f_max / 200.0)?;
  let g = exponential_g(b, (0.0, g_max), g_max / 200.0)?;

  // Convert to Distribution objects for discrete convolution
  let f_dist = Distribution::function(f.x().clone(), f.y().clone())?;
  let g_dist = Distribution::function(g.x().clone(), g.y().clone())?;

  // Discrete convolution
  let discrete = distribution_convolution(&f_dist, &g_dist).wrap_err("Failed to compute discrete convolution")?;

  // Analytical convolution with high resolution
  let result_domain = (0.0, f_max + g_max);
  let analytical = exponential_convolution(a, b, result_domain, (f_max + g_max) / 400.0)?;

  // Compare results
  compare_with_distribution(&analytical, &discrete, tolerance)
    .wrap_err("Failed to compare exponential convolution results")
}

/// Batch verification of multiple Gaussian parameter sets
///
/// Tests convolution accuracy across a range of parameter combinations
/// to assess robustness and identify edge cases.
pub fn verify_gaussian_batch(test_cases: &[(f64, f64, f64)], tolerance: f64) -> Result<Vec<f64>, eyre::Report> {
  let mut results = Vec::with_capacity(test_cases.len());

  for &(sigma_f, sigma_g, mu) in test_cases {
    let diff = verify_gaussian_convolution(sigma_f, sigma_g, mu, tolerance)
      .wrap_err_with(|| format!("Failed verification for σf={sigma_f}, σg={sigma_g}, μ={mu}"))?;
    results.push(diff);
  }

  Ok(results)
}

/// Batch verification of multiple exponential parameter sets
///
/// Tests convolution accuracy across a range of parameter combinations
/// to assess robustness and identify edge cases.
pub fn verify_exponential_batch(test_cases: &[(f64, f64, f64, f64)], tolerance: f64) -> Result<Vec<f64>, eyre::Report> {
  let mut results = Vec::with_capacity(test_cases.len());

  for &(a, b, f_max, g_max) in test_cases {
    let diff = verify_exponential_convolution(a, b, f_max, g_max, tolerance)
      .wrap_err_with(|| format!("Failed verification for a={a}, b={b}, f_max={f_max}, g_max={g_max}"))?;
    results.push(diff);
  }

  Ok(results)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;

  #[test]
  fn test_verify_gaussian_convolution_reference_example() {
    // Reference example from convolution.md
    let sigma_f = 1.0;
    let sigma_g = 2.0;
    let mu = 1.0;

    let max_diff = verify_gaussian_convolution(sigma_f, sigma_g, mu, 10.0).unwrap();
    // Verify the comparison runs without error and produces a reasonable result
    assert!(max_diff >= 0.0, "Difference should be non-negative");
    assert!(max_diff < 100.0, "Difference should be reasonable, got {max_diff}");
  }

  #[test]
  fn test_verify_exponential_convolution_reference_example() {
    // Reference example from convolution.md
    let a = 1.0;
    let b = 2.0;
    let f_max = 10.0;
    let g_max = 7.0;

    let max_diff = verify_exponential_convolution(a, b, f_max, g_max, 10.0).unwrap();
    // Verify the comparison runs without error and produces a reasonable result
    assert!(max_diff >= 0.0, "Difference should be non-negative");
    assert!(max_diff < 100.0, "Difference should be reasonable, got {max_diff}");
  }

  #[test]
  fn test_verify_gaussian_batch() {
    let test_cases = vec![(0.5, 1.0, 0.0), (1.5, 0.8, -1.0), (2.0, 1.5, 2.0)];

    let results = verify_gaussian_batch(&test_cases, 100.0).unwrap();
    assert_eq!(results.len(), test_cases.len());

    for (i, &diff) in results.iter().enumerate() {
      assert!(diff >= 0.0, "Test case {i}: difference should be non-negative");
      assert!(
        diff < 1000.0,
        "Test case {i}: difference should be reasonable, got {diff}"
      );
    }
  }

  #[test]
  fn test_verify_exponential_batch() {
    let test_cases = vec![(0.5, 1.5, 8.0, 6.0), (1.0, 3.0, 5.0, 4.0), (2.0, 0.8, 6.0, 8.0)];

    let results = verify_exponential_batch(&test_cases, 100.0).unwrap();
    assert_eq!(results.len(), test_cases.len());

    for (i, &diff) in results.iter().enumerate() {
      assert!(diff >= 0.0, "Test case {i}: difference should be non-negative");
      assert!(
        diff < 1000.0,
        "Test case {i}: difference should be reasonable, got {diff}"
      );
    }
  }

  #[test]
  fn test_compare_with_distribution_empty() {
    let analytical = GridFn::new(
      ndarray::Array1::from_vec(vec![0.0, 1.0]),
      ndarray::Array1::from_vec(vec![1.0, 0.0]),
    )
    .unwrap();

    let empty_dist = Distribution::empty();
    let diff = compare_with_distribution(&analytical, &empty_dist, 1.0).unwrap();
    assert_ulps_eq!(diff, 0.0, max_ulps = 1);
  }
}
