#![allow(clippy::many_single_char_names)]
use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

use super::test_cases::{ExponentialTestCase, create_exponential_test_cases};

pub struct ExponentialConvInput {
  test_cases: Vec<ExponentialTestCase>,
}

impl ConvInput for ExponentialConvInput {
  type TestCase = ExponentialTestCase;

  fn new() -> Self {
    Self {
      test_cases: create_exponential_test_cases(),
    }
  }

  fn with_test_cases(test_cases: Vec<ExponentialTestCase>) -> Self {
    Self { test_cases }
  }

  fn create_f(&self, &ExponentialTestCase { a, f_domain, dx, .. }: &Self::TestCase) -> Result<GridFn, Report> {
    // Generate exponential function f(x) = a*exp(-ax) for x ≥ 0, 0 otherwise
    GridFn::from_grid(f_domain, dx, |x| if x >= 0.0 { a * (-a * x).exp() } else { 0.0 })
  }

  fn create_g(&self, &ExponentialTestCase { b, g_domain, dx, .. }: &Self::TestCase) -> Result<GridFn, Report> {
    // Generate exponential function g(x) = b*exp(-bx) for x ≥ 0, 0 otherwise
    GridFn::from_grid(g_domain, dx, |x| if x >= 0.0 { b * (-b * x).exp() } else { 0.0 })
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(
    &self,
    &ExponentialTestCase {
      a, b, eval_domain, dx, ..
    }: &Self::TestCase,
    _eval_grid: &Array1<f64>,
  ) -> Result<GridFn, Report> {
    // Analytical convolution of two exponential functions
    //
    // f(x) = a*exp(-ax) for x ≥ 0, 0 otherwise
    // g(x) = b*exp(-bx) for x ≥ 0, 0 otherwise
    //
    // Result: (ab)/(a-b) * (1-exp(-(a-b)x)) * exp(-bx) for a ≠ b
    // Special case: when a = b, result is ab*x*exp(-ax)
    GridFn::from_grid(eval_domain, dx, |x| {
      if x < 0.0 {
        0.0
      } else if (a - b).abs() < 1e-15 {
        a * b * x * (-a * x).exp()
      } else {
        (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp()
      }
    })
  }

  fn test_cases(&self) -> &[Self::TestCase] {
    &self.test_cases
  }

  fn function_type(&self) -> &'static str {
    "exponential"
  }
}
