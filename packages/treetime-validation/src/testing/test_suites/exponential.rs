#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::ConvolutionTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_analytical::validation::cases::exponential_convolution::{
  EXPONENTIAL_CONVOLUTION_CASES, ExponentialConvolutionTestCase,
};
use treetime_analytical::{exponential_convolution_grid, exponential_pdf_grid};

#[derive(Default)]
pub struct ExponentialTestSuite;

impl ConvolutionTestSuite for ExponentialTestSuite {
  type TestCase = ExponentialTestCase;

  fn test_suite_name(&self) -> &'static str {
    "exponential"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(exponential_pdf_grid(test_case.a, grid))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(exponential_pdf_grid(test_case.b, grid))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(exponential_convolution_grid(test_case.a, test_case.b, eval_grid))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    EXPONENTIAL_CONVOLUTION_CASES
      .iter()
      .map(ExponentialTestCase::from)
      .collect()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
  pub a: f64,
  pub b: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl From<&ExponentialConvolutionTestCase> for ExponentialTestCase {
  fn from(case: &ExponentialConvolutionTestCase) -> Self {
    Self {
      name: case.name.to_owned(),
      description: case.description.to_owned(),
      stress_type: case.stress_type.to_owned(),
      analytical_caution: case.analytical_caution.to_owned(),
      slowness: case.slowness,
      a: case.a,
      b: case.b,
      input_grid_domain: case.input_grid_domain,
      input_grid_n_points: case.input_grid_n_points,
    }
  }
}

impl TestCase for ExponentialTestCase {
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

  fn slowness(&self) -> f64 {
    self.slowness
  }

  fn input_grid_domain(&self) -> (f64, f64) {
    self.input_grid_domain
  }

  fn input_grid_n_points(&self) -> usize {
    self.input_grid_n_points
  }
}
