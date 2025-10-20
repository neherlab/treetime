use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

use super::analytical::{gaussian_convolution, gaussian_f, gaussian_g};
use super::test_cases::{GaussianTestCase, create_gaussian_test_cases};

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

  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    gaussian_f(test_case.sigma_f, test_case.f_domain, test_case.dx)
  }

  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    gaussian_g(test_case.sigma_g, test_case.mu, test_case.g_domain, test_case.dx)
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, _eval_grid: &Array1<f64>) -> Result<GridFn, Report> {
    gaussian_convolution(
      test_case.sigma_f,
      test_case.sigma_g,
      test_case.mu,
      test_case.eval_domain,
      test_case.dx,
    )
  }

  fn test_cases(&self) -> &[Self::TestCase] {
    &self.test_cases
  }

  fn function_type(&self) -> &'static str {
    "gaussian"
  }
}
