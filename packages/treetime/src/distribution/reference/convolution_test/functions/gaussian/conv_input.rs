#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::convolution_test::functions::gaussian::test_cases::{
  GaussianTestCase, create_gaussian_test_cases,
};
use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;
use std::f64::consts::PI;

pub struct GaussianConvInput {
  test_cases: Vec<GaussianTestCase>,
}

impl ConvInput for GaussianConvInput {
  type TestCase = GaussianTestCase;

  fn new() -> Self {
    Self {
      test_cases: create_gaussian_test_cases(),
    }
  }

  fn with_test_cases(test_cases: Vec<GaussianTestCase>) -> Self {
    Self { test_cases }
  }

  fn create_f(
    &self,
    &GaussianTestCase {
      sigma_f, f_domain, dx, ..
    }: &Self::TestCase,
  ) -> Result<GridFn, Report> {
    // Generate Gaussian function f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
    GridFn::from_grid(f_domain, dx, |x| {
      (-(0.5 * (x / sigma_f).powi(2))).exp() / (sigma_f * (2.0 * PI).sqrt())
    })
  }

  fn create_g(
    &self,
    &GaussianTestCase {
      sigma_g,
      mu,
      g_domain,
      dx,
      ..
    }: &Self::TestCase,
  ) -> Result<GridFn, Report> {
    // Generate Gaussian function g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
    GridFn::from_grid(g_domain, dx, |x| {
      (-(0.5 * ((x - mu) / sigma_g).powi(2))).exp() / (sigma_g * (2.0 * PI).sqrt())
    })
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(
    &self,
    &GaussianTestCase {
      sigma_f,
      sigma_g,
      mu,
      eval_domain,
      dx,
      ..
    }: &Self::TestCase,
    _eval_grid: &Array1<f64>,
  ) -> Result<GridFn, Report> {
    // Analytical convolution of two Gaussian functions
    //
    // f(x) = (1/√(2π σ_f²)) * exp(-x²/(2σ_f²))
    // g(x) = (1/√(2π σ_g²)) * exp(-(x-μ)²/(2σ_g²))
    //
    // Result: (1/√(2π(σ_f² + σ_g²))) * exp(-(x-μ)²/(2(σ_f² + σ_g²)))
    let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
    GridFn::from_grid(eval_domain, dx, |x| {
      (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()
    })
  }

  fn test_cases(&self) -> &[Self::TestCase] {
    &self.test_cases
  }

  fn function_type(&self) -> &'static str {
    "gaussian"
  }
}
