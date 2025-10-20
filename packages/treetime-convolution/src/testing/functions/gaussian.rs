#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Default)]
pub struct GaussianConvInput;

impl TestSuite for GaussianConvInput {
  type TestCase = GaussianTestCase;

  fn function_type(&self) -> &'static str {
    "gaussian"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let sigma_f = test_case.sigma_f;
    Ok(grid.mapv(|x| (-(0.5 * (x / sigma_f).powi(2))).exp() / (sigma_f * (2.0 * PI).sqrt())))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let sigma_g = test_case.sigma_g;
    let mu = test_case.mu;
    Ok(grid.mapv(|x| (-(0.5 * ((x - mu) / sigma_g).powi(2))).exp() / (sigma_g * (2.0 * PI).sqrt())))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let sigma_f = test_case.sigma_f;
    let sigma_g = test_case.sigma_g;
    let mu = test_case.mu;
    let variance_sum = sigma_f.powi(2) + sigma_g.powi(2);
    Ok(eval_grid.mapv(|x| (-(0.5 * (x - mu).powi(2) / variance_sum)).exp() / (2.0 * PI * variance_sum).sqrt()))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    vec![
      GaussianTestCase {
        name: "python_notebook_case_1".to_owned(),
        description: "Parameters from conv_gauss_py.ipynb: σ_f=1, σ_g=1, μ=1".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 1.0,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 101,
        output_grid_domain: (-9.0, 11.0),
        output_grid_n_points: 101,
      },
      GaussianTestCase {
        name: "python_notebook_case_2".to_owned(),
        description: "Parameters from conv_gauss_py.ipynb: σ_f=2, σ_g=1, μ=1".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 2.0,
        sigma_g: 1.0,
        mu: 1.0,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 101,
        output_grid_domain: (-9.0, 11.0),
        output_grid_n_points: 101,
      },
      GaussianTestCase {
        name: "tight_truncation".to_owned(),
        description: "Tight truncation highlighting sensitivity to insufficient domain coverage.".to_owned(),
        stress_type: "boundary effects, wrap-around".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 2.0,
        mu: 0.0,
        input_grid_domain: (-6.0, 6.0),
        input_grid_n_points: 1201,
        output_grid_domain: (-8.0, 8.0),
        output_grid_n_points: 1601,
      },
      GaussianTestCase {
        name: "coarse_grid".to_owned(),
        description: "Coarse grid coverage emphasizing discretization error with sub-grid shift.".to_owned(),
        stress_type: "aliasing, grid spacing scaling".to_owned(),
        analytical_caution: "none".to_owned(),
        sigma_f: 1.0,
        sigma_g: 1.0,
        mu: 0.5,
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 161,
        output_grid_domain: (-10.0, 10.0),
        output_grid_n_points: 201,
      },
    ]
  }
}

/// Test case parameters for Gaussian convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub sigma_f: f64,
  pub sigma_g: f64,
  pub mu: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
  pub output_grid_domain: (f64, f64),
  pub output_grid_n_points: usize,
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

  fn input_grid_domain(&self) -> (f64, f64) {
    self.input_grid_domain
  }

  fn input_grid_n_points(&self) -> usize {
    self.input_grid_n_points
  }

  fn output_grid_domain(&self) -> (f64, f64) {
    self.output_grid_domain
  }

  fn output_grid_n_points(&self) -> usize {
    self.output_grid_n_points
  }
}
