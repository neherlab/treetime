#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_analytical::validation::cases::gaussian_convolution::{
  GAUSSIAN_CONVOLUTION_CASES, GaussianConvolutionTestCase,
};
use treetime_analytical::{gaussian_convolution_pdf_grid, gaussian_pdf_grid};

#[derive(Default)]
pub struct GaussianTestSuite;

impl TestSuite for GaussianTestSuite {
  type TestCase = GaussianTestCase;

  fn test_suite_name(&self) -> &'static str {
    "gaussian"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_pdf_grid(0.0, test_case.sigma_f, grid))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_pdf_grid(test_case.mu, test_case.sigma_g, grid))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(gaussian_convolution_pdf_grid(
      test_case.sigma_f,
      test_case.sigma_g,
      test_case.mu,
      eval_grid,
    ))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    GAUSSIAN_CONVOLUTION_CASES.iter().map(GaussianTestCase::from).collect()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl From<&GaussianConvolutionTestCase> for GaussianTestCase {
  fn from(case: &GaussianConvolutionTestCase) -> Self {
    Self {
      name: case.name.to_owned(),
      description: case.description.to_owned(),
      stress_type: case.stress_type.to_owned(),
      analytical_caution: case.analytical_caution.to_owned(),
      slowness: case.slowness,
      sigma_f: case.sigma_f,
      sigma_g: case.sigma_g,
      mu: case.mu,
      input_grid_domain: case.input_grid_domain,
      input_grid_n_points: case.input_grid_n_points,
    }
  }
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
