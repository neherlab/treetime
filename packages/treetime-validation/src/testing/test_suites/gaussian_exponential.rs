#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::{TestCase, TestCaseBase};
use crate::testing::test_suites::test_suites::ConvolutionTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use treetime_analytical::validation::cases::gaussian_exponential::{
  GAUSSIAN_EXPONENTIAL_CASES, GaussianExponentialTestCase as AnalyticalCase,
};
use treetime_analytical::{exponential_pdf_grid, gaussian_exponential_convolution_grid};

#[derive(Default)]
pub struct GaussianExponentialTestSuite;

impl ConvolutionTestSuite for GaussianExponentialTestSuite {
  type TestCase = GaussianExponentialTestCase;

  fn test_suite_name(&self) -> &'static str {
    "conv-gaussian-exponential"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(exponential_pdf_grid(test_case.a_f, grid))
  }

  fn create_g(&self, _test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(grid.mapv(|x| (-(0.5 * x.powi(2))).exp() / (2.0 * PI).sqrt()))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_exponential_convolution_grid(test_case.a_f, eval_grid))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    GAUSSIAN_EXPONENTIAL_CASES
      .iter()
      .map(GaussianExponentialTestCase::from)
      .collect()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianExponentialTestCase {
  #[serde(flatten)]
  pub base: TestCaseBase,
  pub a_f: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl From<&AnalyticalCase> for GaussianExponentialTestCase {
  fn from(case: &AnalyticalCase) -> Self {
    Self {
      base: TestCaseBase {
        name: case.name.to_owned(),
        description: case.description.to_owned(),
        stress_type: case.stress_type.to_owned(),
        analytical_caution: case.analytical_caution.to_owned(),
        slowness: case.slowness,
      },
      a_f: case.a_f,
      input_grid_domain: case.input_grid_domain,
      input_grid_n_points: case.input_grid_n_points,
    }
  }
}

impl TestCase for GaussianExponentialTestCase {
  fn base(&self) -> &TestCaseBase {
    &self.base
  }

  fn input_grid_domain(&self) -> (f64, f64) {
    self.input_grid_domain
  }

  fn input_grid_n_points(&self) -> usize {
    self.input_grid_n_points
  }
}
