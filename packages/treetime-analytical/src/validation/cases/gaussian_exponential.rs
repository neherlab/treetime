use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianExponentialTestCase {
  pub name: &'static str,
  pub description: &'static str,
  pub stress_type: &'static str,
  pub analytical_caution: &'static str,
  pub slowness: f64,
  pub a_f: f64,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

pub const GAUSSIAN_EXPONENTIAL_CASES: &[GaussianExponentialTestCase] = &[
  GaussianExponentialTestCase {
    name: "python_notebook_case",
    description: "Parameters from conv_gauss_exp.ipynb: a_f=0.5",
    stress_type: "reference implementation validation",
    analytical_caution: "none",
    slowness: 0.0,
    a_f: 0.5,
    input_grid_domain: (-5.0, 40.0),
    input_grid_n_points: 101,
  },
  GaussianExponentialTestCase {
    name: "python_notebook_case_fine",
    description: "Parameters from conv_gauss_exp.ipynb: a_f=0.5",
    stress_type: "reference implementation validation",
    analytical_caution: "none",
    slowness: 0.0,
    a_f: 0.5,
    input_grid_domain: (-5.0, 40.0),
    input_grid_n_points: 501,
  },
];
