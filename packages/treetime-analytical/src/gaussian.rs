use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use treetime_ops::ScaledArray;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
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
///
/// # Preconditions
///
/// All `sigma` values must be positive (non-zero). Zero sigma represents a Dirac delta
/// which cannot be represented as a Gaussian product.
pub fn gaussian_product_params(params: &[GaussianParams]) -> (f64, f64, f64) {
  debug_assert!(
    params.iter().all(|p| p.sigma > 0.0),
    "gaussian_product_params: all sigma values must be positive"
  );

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
/// Returns `ScaledArray { normalized, log_scale }` where the full result is:
/// `normalized * exp(log_scale)`.
pub fn gaussian_product(params: &[GaussianParams], grid: &Array1<f64>) -> ScaledArray {
  if params.is_empty() {
    return ScaledArray::new(Array1::ones(grid.len()), 0.0);
  }

  let (mu_star, sigma_star, log_scale) = gaussian_product_params(params);
  let normalized = grid.mapv(|x| (-(0.5 * ((x - mu_star) / sigma_star).powi(2))).exp());
  ScaledArray::new(normalized, log_scale)
}

/// Evaluate a single Gaussian on grid.
pub fn gaussian_evaluate(params: &GaussianParams, grid: &Array1<f64>) -> Array1<f64> {
  grid.mapv(|x| params.amplitude * (-(0.5 * ((x - params.mu) / params.sigma).powi(2))).exp())
}

/// Normalized Gaussian PDF.
///
/// f(x) = exp(-0.5 * ((x - mu) / sigma)^2) / (sigma * sqrt(2 * pi))
///
/// # Preconditions
///
/// `sigma` must be positive. Zero sigma represents a Dirac delta which has no finite PDF value.
pub fn gaussian_pdf(mu: f64, sigma: f64, x: f64) -> f64 {
  debug_assert!(sigma > 0.0, "gaussian_pdf: sigma must be positive");
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
