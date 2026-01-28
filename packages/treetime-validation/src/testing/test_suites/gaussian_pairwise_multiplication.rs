#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::MultiplicationTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_analytical::validation::cases::gaussian_pairwise_multiplication::{
  GAUSSIAN_PAIRWISE_MULTIPLICATION_CASES, GaussianPairwiseMultiplicationTestCase as AnalyticalCase,
};
use treetime_analytical::{GaussianParams, gaussian_product};
use treetime_ops::ScaledArray;

#[derive(Default)]
pub struct GaussianPairwiseMultiplicationTestSuite;

impl MultiplicationTestSuite for GaussianPairwiseMultiplicationTestSuite {
  type TestCase = GaussianPairwiseMultiplicationTestCase;

  fn test_suite_name(&self) -> &'static str {
    "mult-gaussian-pairwise"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let mu = test_case.mu_f;
    let sigma = test_case.sigma_f;
    let amplitude = test_case.amplitude_f;
    Ok(grid.mapv(|x| amplitude * (-(0.5 * ((x - mu) / sigma).powi(2))).exp()))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let mu = test_case.mu_g;
    let sigma = test_case.sigma_g;
    let amplitude = test_case.amplitude_g;
    Ok(grid.mapv(|x| amplitude * (-(0.5 * ((x - mu) / sigma).powi(2))).exp()))
  }

  fn analytical_multiplication(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<ScaledArray, Report> {
    let params = [
      GaussianParams {
        mu: test_case.mu_f,
        sigma: test_case.sigma_f,
        amplitude: test_case.amplitude_f,
      },
      GaussianParams {
        mu: test_case.mu_g,
        sigma: test_case.sigma_g,
        amplitude: test_case.amplitude_g,
      },
    ];

    Ok(gaussian_product(&params, grid))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    GAUSSIAN_PAIRWISE_MULTIPLICATION_CASES
      .iter()
      .map(GaussianPairwiseMultiplicationTestCase::from)
      .collect()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianPairwiseMultiplicationTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
  pub mu_f: f64,
  pub sigma_f: f64,
  pub amplitude_f: f64,
  pub mu_g: f64,
  pub sigma_g: f64,
  pub amplitude_g: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl From<&AnalyticalCase> for GaussianPairwiseMultiplicationTestCase {
  fn from(case: &AnalyticalCase) -> Self {
    Self {
      name: case.name.to_owned(),
      description: case.description.to_owned(),
      stress_type: case.stress_type.to_owned(),
      analytical_caution: case.analytical_caution.to_owned(),
      slowness: case.slowness,
      mu_f: case.mu_f,
      sigma_f: case.sigma_f,
      amplitude_f: case.amplitude_f,
      mu_g: case.mu_g,
      sigma_g: case.sigma_g,
      amplitude_g: case.amplitude_g,
      input_grid_domain: case.input_grid_domain,
      input_grid_n_points: case.input_grid_n_points,
    }
  }
}

impl TestCase for GaussianPairwiseMultiplicationTestCase {
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
