use crate::distribution::reference::convolution_test::traits::ConvInput;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

use super::analytical::{exponential_convolution, exponential_f, exponential_g};
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

  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    exponential_f(test_case.a, test_case.f_domain, test_case.dx)
  }

  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report> {
    exponential_g(test_case.b, test_case.g_domain, test_case.dx)
  }

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64) {
    test_case.eval_domain
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, _eval_grid: &Array1<f64>) -> Result<GridFn, Report> {
    exponential_convolution(test_case.a, test_case.b, test_case.eval_domain, test_case.dx)
  }

  fn test_cases(&self) -> &[Self::TestCase] {
    &self.test_cases
  }

  fn function_type(&self) -> &'static str {
    "exponential"
  }
}
