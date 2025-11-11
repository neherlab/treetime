#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use statrs::function::erf::erfc;
use std::f64::consts::PI;

#[derive(Default)]
pub struct GaussianExponentialTestSuite;

impl TestSuite for GaussianExponentialTestSuite {
  type TestCase = GaussianExponentialTestCase;

  fn test_suite_name(&self) -> &'static str {
    "gaussian-exponential"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let a_f = test_case.a_f;
    Ok(grid.mapv(|x| if x >= 0.0 { a_f * (-a_f * x).exp() } else { 0.0 }))
  }

  fn create_g(&self, _test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    Ok(grid.mapv(|x| (-(0.5 * x.powi(2))).exp() / (2.0 * PI).sqrt()))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let a_f = test_case.a_f;
    Ok(eval_grid.mapv(|x| 0.5 * a_f * (-x * a_f + 0.5 * a_f.powi(2)).exp() * erfc((a_f - x) / 2_f64.sqrt())))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    vec![
      GaussianExponentialTestCase {
        name: "python_notebook_case".to_owned(),
        description: "Parameters from conv_gauss_exp.ipynb: a_f=0.5".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a_f: 0.5,
        input_grid_domain: (-5.0, 40.0),
        input_grid_n_points: 101,
      },
      GaussianExponentialTestCase {
        name: "python_notebook_case_fine".to_owned(),
        description: "Parameters from conv_gauss_exp.ipynb: a_f=0.5".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        a_f: 0.5,
        input_grid_domain: (-5.0, 40.0),
        input_grid_n_points: 501,
      },
    ]
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianExponentialTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
  pub a_f: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl TestCase for GaussianExponentialTestCase {
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
