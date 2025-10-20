use crate::testing::framework::test_case::TestCase;
use eyre::Report;
use ndarray::Array1;

pub trait TestSuite: Send + Sync {
  type TestCase: TestCase;

  fn function_type(&self) -> &'static str;

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn create_test_cases(&self) -> Vec<Self::TestCase>;
}
