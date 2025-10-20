#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Default)]
pub struct GaussianConvInput;

impl GaussianConvInput {
  fn create_all_test_cases() -> Vec<GaussianTestCase> {
    vec![
      GaussianTestCase {
        name: "tight_truncation".to_owned(),
        description: "Tight truncation highlighting sensitivity to insufficient domain coverage.".to_owned(),
        stress_type: "boundary effects, wrap-around".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 2.0,
        mu: 0.0,
        f_domain: (-3.0, 3.0),
        g_domain: (-6.0, 6.0),
        eval_domain: (-4.0, 4.0),
        dx: 0.01,
      },
      GaussianTestCase {
        name: "coarse_grid".to_owned(),
        description: "Coarse grid coverage emphasizing discretization error with sub-grid shift.".to_owned(),
        stress_type: "aliasing, grid spacing scaling".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 0.5,
        f_domain: (-8.0, 8.0),
        g_domain: (-7.5, 8.5),
        eval_domain: (-6.0, 6.0),
        dx: 0.1,
      },
    ]
  }
}

impl ConvInput for GaussianConvInput {
  type TestCase = GaussianTestCase;

  fn function_type(&self) -> &'static str {
    "gaussian"
  }

  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    let GaussianTestCase {
      sigma_f, f_domain, dx, ..
    } = test_case;
    // Generate Gaussian function f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
    GridFn::from_grid(*f_domain, *dx, |x| {
      (-(0.5 * (x / sigma_f).powi(2))).exp() / (sigma_f * (2.0 * PI).sqrt())
    })
  }

  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    let GaussianTestCase {
      sigma_g,
      mu,
      g_domain,
      dx,
      ..
    } = test_case;
    // Generate Gaussian function g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
    GridFn::from_grid(*g_domain, *dx, |x| {
      (-(0.5 * ((x - mu) / sigma_g).powi(2))).exp() / (sigma_g * (2.0 * PI).sqrt())
    })
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, _eval_grid: &Array1<f64>) -> Result<GridFn, Report> {
    let GaussianTestCase {
      sigma_f,
      sigma_g,
      mu,
      eval_domain,
      dx,
      ..
    } = test_case;
    // Analytical convolution of two Gaussian functions
    //
    // f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
    // g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
    //
    // Result: (1/√(2π(σ_f² + σ_g²))) * exp(-(x-μ)²/(2(σ_f² + σ_g²)))
    let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
    GridFn::from_grid(*eval_domain, *dx, |x| {
      (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()
    })
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    Self::create_all_test_cases()
  }
}

/// Test case parameters for Gaussian convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub f_domain: (f64, f64),
  pub g_domain: (f64, f64),
  pub eval_domain: (f64, f64),
  pub dx: f64,
}

impl TestCase for GaussianTestCase {
  fn name(&self) -> &str {
    &self.name
  }

  fn description(&self) -> &str {
    &self.description
  }

  fn stress_type(&self) -> &str {
    &self.stress_type
  }

  fn analytical_caution(&self) -> &str {
    &self.analytical_caution
  }

  fn dx(&self) -> f64 {
    self.dx
  }
}
