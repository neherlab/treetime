use super::*;
use std::f64::consts::PI;

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

  let result = gaussian_product(&params, &grid);

  assert_relative_eq!(result.normalized[2], 1.0, epsilon = 1e-10);
  assert!(result.normalized[0] < result.normalized[2]);
  assert!(result.normalized[4] < result.normalized[2]);
  assert_relative_eq!(result.normalized[0], result.normalized[4], epsilon = 1e-10);
  assert_relative_eq!(result.normalized[1], result.normalized[3], epsilon = 1e-10);
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
  let expected = 1.0 / (sigma * (2.0_f64 * PI).sqrt());
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
  let expected = 1.0 / (expected_sigma * (2.0_f64 * PI).sqrt());
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
