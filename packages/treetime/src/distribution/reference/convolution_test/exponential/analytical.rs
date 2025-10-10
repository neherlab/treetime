#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

/// Generate exponential function f(x) = a*exp(-ax) for x ≥ 0, 0 otherwise
pub fn exponential_f(a: f64, domain: (f64, f64), dx: f64) -> Result<GridFn, Report> {
  GridFn::from_grid(domain, dx, |x| if x >= 0.0 { a * (-a * x).exp() } else { 0.0 })
}

/// Generate exponential function g(x) = b*exp(-bx) for x ≥ 0, 0 otherwise
pub fn exponential_g(b: f64, domain: (f64, f64), dx: f64) -> Result<GridFn, Report> {
  GridFn::from_grid(domain, dx, |x| if x >= 0.0 { b * (-b * x).exp() } else { 0.0 })
}

/// Generate exponential function f(x) = a*exp(-ax) on shared grid
pub fn exponential_f_on_grid(a: f64, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  GridFn::from_shared_grid(x_grid, |x| if x >= 0.0 { a * (-a * x).exp() } else { 0.0 })
}

/// Generate exponential function g(x) = b*exp(-bx) on shared grid
pub fn exponential_g_on_grid(b: f64, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  GridFn::from_shared_grid(x_grid, |x| if x >= 0.0 { b * (-b * x).exp() } else { 0.0 })
}

/// Generate analytical convolution on shared grid
pub fn exponential_convolution_on_grid(a: f64, b: f64, x_grid: &Array1<f64>) -> Result<GridFn, Report> {
  GridFn::from_shared_grid(x_grid, |x| {
    if x < 0.0 {
      0.0
    } else if (a - b).abs() < 1e-15 {
      a * b * x * (-a * x).exp()
    } else {
      (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp()
    }
  })
}

/// Analytical convolution of two exponential functions
///
/// f(x) = a*exp(-ax) for x ≥ 0, 0 otherwise
/// g(x) = b*exp(-bx) for x ≥ 0, 0 otherwise
///
/// Result: (ab)/(a-b) * (1-exp(-(a-b)x)) * exp(-bx) for a ≠ b
/// Special case: when a = b, result is ab*x*exp(-ax)
pub fn exponential_convolution(a: f64, b: f64, domain: (f64, f64), dx: f64) -> Result<GridFn, Report> {
  GridFn::from_grid(domain, dx, |x| {
    if x < 0.0 {
      0.0
    } else if (a - b).abs() < 1e-15 {
      a * b * x * (-a * x).exp()
    } else {
      (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp()
    }
  })
}
