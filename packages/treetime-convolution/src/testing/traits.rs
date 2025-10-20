use crate::grid_fn::GridFn;
use crate::testing::framework::test_case::TestCase;
use eyre::Report;
use ndarray::Array1;

pub trait TestSuite: Send + Sync {
  type TestCase: TestCase;

  fn function_type(&self) -> &'static str;

  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64);

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<GridFn, Report>;

  fn create_test_cases(&self) -> Vec<Self::TestCase>;
}

pub trait Algo: Send + Sync {
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;

  fn name(&self) -> &'static str;
}
