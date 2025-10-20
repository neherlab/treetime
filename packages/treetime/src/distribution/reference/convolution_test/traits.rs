use crate::distribution::reference::convolution_test::framework::test_case::TestCase;
use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

pub trait ConvInput: Send + Sync {
  type TestCase: TestCase;

  fn function_type(&self) -> &'static str;

  fn create_f(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  fn create_g(&self, test_case: &Self::TestCase) -> Result<GridFn, Report>;

  fn eval_domain(&self, test_case: &Self::TestCase) -> (f64, f64);

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<GridFn, Report>;

  fn create_test_cases(&self) -> Vec<Self::TestCase>;
}

pub trait ConvAlgo: Send + Sync {
  fn convolve(&self, f: &GridFn, g: &GridFn, x_grid: &Array1<f64>) -> Result<GridFn, Report>;

  fn name(&self) -> &'static str;
}
