#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::TestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Default)]
pub struct ExponentialConvInput;

impl TestSuite for ExponentialConvInput {
  type TestCase = ExponentialTestCase;

  fn function_type(&self) -> &'static str {
    "exponential"
  }

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let a = test_case.a;
    Ok(grid.mapv(|x| if x >= 0.0 { a * (-a * x).exp() } else { 0.0 }))
  }

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let b = test_case.b;
    Ok(grid.mapv(|x| if x >= 0.0 { b * (-b * x).exp() } else { 0.0 }))
  }

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let a = test_case.a;
    let b = test_case.b;
    Ok(eval_grid.mapv(|x| {
      if x < 0.0 {
        0.0
      } else if (a - b).abs() < 1e-15 {
        a * b * x * (-a * x).exp()
      } else {
        (a * b) / (a - b) * (1.0 - (-(a - b) * x).exp()) * (-b * x).exp()
      }
    }))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    vec![
      ExponentialTestCase {
        name: "python_notebook_case".to_owned(),
        description: "Parameters from conv_exp_py.ipynb: a=1, b=2".to_owned(),
        stress_type: "reference implementation validation".to_owned(),
        analytical_caution: "none".to_owned(),
        a: 1.0,
        b: 2.0,
        input_grid_domain: (-1.0, 10.0),
        input_grid_n_points: 1101,
        output_grid_domain: (-1.0, 10.0),
        output_grid_n_points: 1101,
      },
      ExponentialTestCase {
        name: "moderate_coarse_grid".to_owned(),
        description: "Moderate rates on a coarse grid. Discretization error on coarse grids; step size sensitivity."
          .to_owned(),
        stress_type: "aliasing/accuracy vs Δx, Δx scaling factor".to_owned(),
        analytical_caution: "none".to_owned(),
        a: 1.0,
        b: 0.8,
        input_grid_domain: (0.0, 20.0),
        input_grid_n_points: 201,
        output_grid_domain: (0.0, 30.0),
        output_grid_n_points: 301,
      },
      ExponentialTestCase {
        name: "tight_truncation".to_owned(),
        description:
          "Tight truncation to expose intentional errors. Sensitivity to insufficient support; boundary effects."
            .to_owned(),
        stress_type: "wrap-around artifacts if padding insufficient".to_owned(),
        analytical_caution: "none".to_owned(),
        a: 1.0,
        b: 2.0,
        input_grid_domain: (0.0, 5.0),
        input_grid_n_points: 501,
        output_grid_domain: (0.0, 6.0),
        output_grid_n_points: 601,
      },
    ]
  }
}

/// Test case parameters for Exponential convolution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExponentialTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub a: f64,
  pub b: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
  pub output_grid_domain: (f64, f64),
  pub output_grid_n_points: usize,
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
