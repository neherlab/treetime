#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use std::f64::consts::PI;

/// Generate Gaussian function f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
pub fn gaussian_f(sigma_f: f64, domain: (f64, f64), dx: f64) -> Result<GridFn, Report> {
  GridFn::from_grid(domain, dx, |x| {
    (-(0.5 * (x / sigma_f).powi(2))).exp() / (sigma_f * (2.0 * PI).sqrt())
  })
}

/// Generate Gaussian function g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
pub fn gaussian_g(sigma_g: f64, mu: f64, domain: (f64, f64), dx: f64) -> Result<GridFn, Report> {
  GridFn::from_grid(domain, dx, |x| {
    (-(0.5 * ((x - mu) / sigma_g).powi(2))).exp() / (sigma_g * (2.0 * PI).sqrt())
  })
}

/// Analytical convolution of two Gaussian functions
///
/// f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
/// g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
///
/// Result: (1/√(2π(σ_f² + σ_g²))) * exp(-(x-μ)²/(2(σ_f² + σ_g²)))
pub fn gaussian_convolution(
  sigma_f: f64,
  sigma_g: f64,
  mu: f64,
  domain: (f64, f64),
  dx: f64,
) -> Result<GridFn, Report> {
  let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
  GridFn::from_grid(domain, dx, |x| {
    (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()
  })
}
