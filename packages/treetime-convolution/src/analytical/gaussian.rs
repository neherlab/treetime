use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GaussianParams {
  pub mu: f64,
  pub sigma: f64,
  pub amplitude: f64,
}

/// Analytical parameters for the product of N Gaussians.
///
/// Product of N Gaussians: G_i(x) = A_i * exp(-0.5 * ((x - mu_i) / sigma_i)^2)
///
/// Result is also Gaussian with:
/// - precision_sum = sum(1 / sigma_i^2)
/// - sigma_star = sqrt(1 / precision_sum)
/// - mu_star = sigma_star^2 * sum(mu_i / sigma_i^2)
/// - quadratic_term = sum((mu_i - mu_star)^2 / sigma_i^2)
/// - log_scale = -0.5 * quadratic_term + sum(ln(A_i))
///
/// Returns (mu_star, sigma_star, log_scale).
pub fn gaussian_product_params(params: &[GaussianParams]) -> (f64, f64, f64) {
  if params.is_empty() {
    return (0.0, f64::INFINITY, 0.0);
  }

  let precision_sum: f64 = params.iter().map(|p| 1.0 / p.sigma.powi(2)).sum();
  let sigma_star = (1.0 / precision_sum).sqrt();
  let mu_star = sigma_star.powi(2) * params.iter().map(|p| p.mu / p.sigma.powi(2)).sum::<f64>();

  let quadratic_term: f64 = params.iter().map(|p| (p.mu - mu_star).powi(2) / p.sigma.powi(2)).sum();
  let log_amplitude_sum: f64 = params.iter().map(|p| p.amplitude.ln()).sum();
  let log_scale = -0.5 * quadratic_term + log_amplitude_sum;

  (mu_star, sigma_star, log_scale)
}

/// Evaluate Gaussian product on grid.
///
/// Returns (normalized_shape, log_scale) where the full result is:
/// `normalized_shape * exp(log_scale)`.
pub fn gaussian_product(params: &[GaussianParams], grid: &Array1<f64>) -> (Array1<f64>, f64) {
  if params.is_empty() {
    return (Array1::ones(grid.len()), 0.0);
  }

  let (mu_star, sigma_star, log_scale) = gaussian_product_params(params);
  let values = grid.mapv(|x| (-(0.5 * ((x - mu_star) / sigma_star).powi(2))).exp());
  (values, log_scale)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_relative_eq;
  use ndarray::array;

  #[test]
  fn test_gaussian_product_params_single() {
    let params = [GaussianParams {
      mu: 1.0,
      sigma: 2.0,
      amplitude: 3.0,
    }];

    let (mu_star, sigma_star, log_scale) = gaussian_product_params(&params);

    assert_relative_eq!(mu_star, 1.0, epsilon = 1e-10);
    assert_relative_eq!(sigma_star, 2.0, epsilon = 1e-10);
    assert_relative_eq!(log_scale, 3.0_f64.ln(), epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_product_params_two_identical() {
    let params = [
      GaussianParams {
        mu: 0.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
      GaussianParams {
        mu: 0.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
    ];

    let (mu_star, sigma_star, log_scale) = gaussian_product_params(&params);

    assert_relative_eq!(mu_star, 0.0, epsilon = 1e-10);
    assert_relative_eq!(sigma_star, 1.0 / 2.0_f64.sqrt(), epsilon = 1e-10);
    assert_relative_eq!(log_scale, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_product_params_two_shifted() {
    let params = [
      GaussianParams {
        mu: -1.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
      GaussianParams {
        mu: 1.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
    ];

    let (mu_star, sigma_star, _log_scale) = gaussian_product_params(&params);

    assert_relative_eq!(mu_star, 0.0, epsilon = 1e-10);
    assert_relative_eq!(sigma_star, 1.0 / 2.0_f64.sqrt(), epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_product_params_empty() {
    let params: [GaussianParams; 0] = [];
    let (mu_star, sigma_star, log_scale) = gaussian_product_params(&params);

    assert_relative_eq!(mu_star, 0.0, epsilon = 1e-10);
    assert!(sigma_star.is_infinite());
    assert_relative_eq!(log_scale, 0.0, epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_product_shape() {
    let params = [
      GaussianParams {
        mu: 0.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
      GaussianParams {
        mu: 0.0,
        sigma: 1.0,
        amplitude: 1.0,
      },
    ];
    let grid = array![-2.0, -1.0, 0.0, 1.0, 2.0];

    let (shape, _log_scale) = gaussian_product(&params, &grid);

    assert_relative_eq!(shape[2], 1.0, epsilon = 1e-10);
    assert!(shape[0] < shape[2]);
    assert!(shape[4] < shape[2]);
    assert_relative_eq!(shape[0], shape[4], epsilon = 1e-10);
    assert_relative_eq!(shape[1], shape[3], epsilon = 1e-10);
  }
}
