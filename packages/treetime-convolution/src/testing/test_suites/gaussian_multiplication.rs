#![allow(clippy::many_single_char_names)]
use crate::analytical::gaussian::{GaussianParams, gaussian_product};
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::MultiplicationTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Default)]
pub struct GaussianMultiplicationTestSuite;

impl MultiplicationTestSuite for GaussianMultiplicationTestSuite {
  type TestCase = GaussianMultiplicationTestCase;

  fn test_suite_name(&self) -> &'static str {
    "gaussian-multiplication"
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

  fn analytical_multiplication(
    &self,
    test_case: &Self::TestCase,
    grid: &Array1<f64>,
  ) -> Result<(Array1<f64>, f64), Report> {
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
    vec![
      GaussianMultiplicationTestCase {
        name: "identical".to_owned(),
        description: "Two identical Gaussians at origin with unit amplitude.".to_owned(),
        stress_type: "baseline".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        mu_f: 0.0,
        sigma_f: 1.0,
        amplitude_f: 1.0,
        mu_g: 0.0,
        sigma_g: 1.0,
        amplitude_g: 1.0,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianMultiplicationTestCase {
        name: "shifted_equal".to_owned(),
        description: "Two Gaussians with equal widths but shifted means; product peak at midpoint.".to_owned(),
        stress_type: "translation".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        mu_f: -1.0,
        sigma_f: 1.0,
        amplitude_f: 1.0,
        mu_g: 1.0,
        sigma_g: 1.0,
        amplitude_g: 1.0,
        input_grid_domain: (-6.0, 6.0),
        input_grid_n_points: 241,
      },
      GaussianMultiplicationTestCase {
        name: "different_widths".to_owned(),
        description: "Two Gaussians with different widths; narrower Gaussian dominates result.".to_owned(),
        stress_type: "width mixing".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        mu_f: 0.0,
        sigma_f: 1.0,
        amplitude_f: 1.0,
        mu_g: 0.0,
        sigma_g: 2.0,
        amplitude_g: 1.0,
        input_grid_domain: (-10.0, 10.0),
        input_grid_n_points: 401,
      },
      GaussianMultiplicationTestCase {
        name: "narrow_wide".to_owned(),
        description: "Extreme width ratio with very narrow and very wide Gaussian.".to_owned(),
        stress_type: "extreme ratio, resolution".to_owned(),
        analytical_caution: "grid resolution may be insufficient for narrow Gaussian".to_owned(),
        slowness: 0.3,
        mu_f: 0.0,
        sigma_f: 0.1,
        amplitude_f: 1.0,
        mu_g: 0.0,
        sigma_g: 10.0,
        amplitude_g: 1.0,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 1001,
      },
      GaussianMultiplicationTestCase {
        name: "shifted_different_widths".to_owned(),
        description: "Shifted means with different widths; tests combined effects.".to_owned(),
        stress_type: "combined".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.1,
        mu_f: -2.0,
        sigma_f: 1.5,
        amplitude_f: 1.0,
        mu_g: 1.0,
        sigma_g: 0.8,
        amplitude_g: 1.0,
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 321,
      },
      GaussianMultiplicationTestCase {
        name: "small_amplitudes".to_owned(),
        description: "Small amplitudes testing log-scale handling.".to_owned(),
        stress_type: "underflow".to_owned(),
        analytical_caution: "numerical underflow possible with plain multiplication".to_owned(),
        slowness: 0.1,
        mu_f: 0.0,
        sigma_f: 1.0,
        amplitude_f: 0.001,
        mu_g: 0.0,
        sigma_g: 1.0,
        amplitude_g: 0.001,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianMultiplicationTestCase {
        name: "large_shift".to_owned(),
        description: "Large shift between means; product has low amplitude.".to_owned(),
        stress_type: "near-zero product".to_owned(),
        analytical_caution: "product may underflow".to_owned(),
        slowness: 0.2,
        mu_f: -5.0,
        sigma_f: 1.0,
        amplitude_f: 1.0,
        mu_g: 5.0,
        sigma_g: 1.0,
        amplitude_g: 1.0,
        input_grid_domain: (-10.0, 10.0),
        input_grid_n_points: 401,
      },
      GaussianMultiplicationTestCase {
        name: "asymmetric_shift".to_owned(),
        description: "Asymmetric shift with different widths.".to_owned(),
        stress_type: "asymmetry".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.15,
        mu_f: 0.0,
        sigma_f: 2.0,
        amplitude_f: 1.0,
        mu_g: 3.0,
        sigma_g: 1.0,
        amplitude_g: 1.0,
        input_grid_domain: (-8.0, 12.0),
        input_grid_n_points: 401,
      },
      GaussianMultiplicationTestCase {
        name: "fine_grid".to_owned(),
        description: "High-resolution grid for reference accuracy.".to_owned(),
        stress_type: "performance".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.5,
        mu_f: 0.0,
        sigma_f: 1.0,
        amplitude_f: 1.0,
        mu_g: 0.5,
        sigma_g: 1.5,
        amplitude_g: 1.0,
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 3201,
      },
      GaussianMultiplicationTestCase {
        name: "mixed_amplitudes".to_owned(),
        description: "Different amplitudes testing scale handling.".to_owned(),
        stress_type: "scale mixing".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.1,
        mu_f: 0.0,
        sigma_f: 1.0,
        amplitude_f: 10.0,
        mu_g: 0.0,
        sigma_g: 1.0,
        amplitude_g: 0.5,
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
    ]
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianMultiplicationTestCase {
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

impl TestCase for GaussianMultiplicationTestCase {
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

