use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

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

/// Evaluate a single Gaussian on grid.
pub fn gaussian_evaluate(params: &GaussianParams, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| params.amplitude * (-(0.5 * ((x - params.mu) / params.sigma).powi(2))).exp())
}

/// Normalized Gaussian PDF.
///
/// f(x) = exp(-0.5 * ((x - mu) / sigma)^2) / (sigma * sqrt(2 * pi))
pub fn gaussian_pdf(mu: f64, sigma: f64, x: f64) -> f64 {
  (-(0.5 * ((x - mu) / sigma).powi(2))).exp() / (sigma * (2.0 * PI).sqrt())
}

/// Evaluate normalized Gaussian PDF on grid.
pub fn gaussian_pdf_grid(mu: f64, sigma: f64, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| gaussian_pdf(mu, sigma, x))
}

/// Analytical convolution of two normalized Gaussian PDFs.
///
/// Given:
/// - f(x) = N(0, sigma_f) centered at origin
/// - g(x) = N(mu, sigma_g) centered at mu
///
/// The convolution is N(mu, sqrt(sigma_f^2 + sigma_g^2)).
pub fn gaussian_convolution_pdf(sigma_f: f64, sigma_g: f64, mu: f64, x: f64) -> f64 {
  let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
  (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()
}

/// Evaluate normalized Gaussian convolution PDF on grid.
pub fn gaussian_convolution_pdf_grid(sigma_f: f64, sigma_g: f64, mu: f64, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| gaussian_convolution_pdf(sigma_f, sigma_g, mu, x))
}

/// Analytical convolution of two Gaussians.
///
/// The convolution of two Gaussians is also a Gaussian with:
/// - mu = mu1 + mu2
/// - sigma = sqrt(sigma1^2 + sigma2^2)
/// - amplitude = amplitude1 * amplitude2 * sqrt(2*pi) * sigma1 * sigma2 / sigma
pub fn gaussian_convolution(a: &GaussianParams, b: &GaussianParams, grid: &Array1<f64>) -> Array1<f64> {
  let sigma_conv = a.sigma.hypot(b.sigma);
  let mu_conv = a.mu + b.mu;

  let normalization = a.amplitude * b.amplitude * std::f64::consts::TAU.sqrt() * a.sigma * b.sigma / sigma_conv;

  grid.mapv(|x| normalization * (-(0.5 * ((x - mu_conv) / sigma_conv).powi(2))).exp())
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

  #[test]
  fn test_gaussian_evaluate() {
    let params = GaussianParams {
      mu: 0.0,
      sigma: 1.0,
      amplitude: 1.0,
    };
    let grid = array![-1.0, 0.0, 1.0];
    let values = gaussian_evaluate(&params, &grid);

    assert_relative_eq!(values[1], 1.0, epsilon = 1e-10);
    assert_relative_eq!(values[0], values[2], epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_convolution_same_width() {
    let a = GaussianParams {
      mu: 0.0,
      sigma: 1.0,
      amplitude: 1.0,
    };
    let b = GaussianParams {
      mu: 0.0,
      sigma: 1.0,
      amplitude: 1.0,
    };
    let grid = array![-3.0, 0.0, 3.0];
    let conv = gaussian_convolution(&a, &b, &grid);

    // Peak should be at center (mu=0)
    assert!(conv[1] > conv[0]);
    assert!(conv[1] > conv[2]);
  }

  #[test]
  fn test_gaussian_pdf_at_mean() {
    let mu = 0.0;
    let sigma = 1.0;
    let result = gaussian_pdf(mu, sigma, mu);
    let expected = 1.0 / (sigma * (2.0 * PI).sqrt());
    assert_relative_eq!(result, expected, epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_pdf_symmetry() {
    let mu = 1.0;
    let sigma = 2.0;
    let delta = 0.5;
    let left = gaussian_pdf(mu, sigma, mu - delta);
    let right = gaussian_pdf(mu, sigma, mu + delta);
    assert_relative_eq!(left, right, epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_pdf_grid() {
    let mu = 0.0;
    let sigma = 1.0;
    let grid = array![-1.0, 0.0, 1.0];
    let result = gaussian_pdf_grid(mu, sigma, &grid);

    assert_relative_eq!(result[0], result[2], epsilon = 1e-10);
    assert!(result[1] > result[0]);
  }

  #[test]
  fn test_gaussian_convolution_pdf_centered() {
    let sigma_f = 1.0;
    let sigma_g = 1.0;
    let mu = 0.0;
    let result = gaussian_convolution_pdf(sigma_f, sigma_g, mu, 0.0);
    let expected_sigma = sigma_f.hypot(sigma_g);
    let expected = 1.0 / (expected_sigma * (2.0 * PI).sqrt());
    assert_relative_eq!(result, expected, epsilon = 1e-10);
  }

  #[test]
  fn test_gaussian_convolution_pdf_shifted() {
    let sigma_f = 1.0;
    let sigma_g = 2.0;
    let mu = 1.0;
    let result_at_mu = gaussian_convolution_pdf(sigma_f, sigma_g, mu, mu);
    let result_away = gaussian_convolution_pdf(sigma_f, sigma_g, mu, mu + 3.0);
    assert!(result_at_mu > result_away);
  }
}
