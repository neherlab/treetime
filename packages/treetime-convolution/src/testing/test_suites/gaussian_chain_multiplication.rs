#![allow(clippy::many_single_char_names)]
use crate::analytical::gaussian::{GaussianParams, gaussian_product};
use crate::testing::framework::test_case::TestCase;
use crate::testing::test_suites::test_suites::ChainMultiplicationTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

#[derive(Default)]
pub struct GaussianChainMultiplicationTestSuite;

impl ChainMultiplicationTestSuite for GaussianChainMultiplicationTestSuite {
  type TestCase = GaussianChainTestCase;

  fn test_suite_name(&self) -> &'static str {
    "gaussian-chain-multiplication"
  }

  fn create_factors(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Vec<Array1<f64>>, Report> {
    let factors = test_case
      .factors
      .iter()
      .map(|p| grid.mapv(|x| p.amplitude * (-(0.5 * ((x - p.mu) / p.sigma).powi(2))).exp()))
      .collect();
    Ok(factors)
  }

  fn analytical_chain_multiplication(
    &self,
    test_case: &Self::TestCase,
    grid: &Array1<f64>,
  ) -> Result<(Array1<f64>, f64), Report> {
    Ok(gaussian_product(&test_case.factors, grid))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    vec![
      GaussianChainTestCase {
        name: "two_identical".to_owned(),
        description: "Two identical Gaussians (baseline for chain).".to_owned(),
        stress_type: "baseline".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        factors: vec![
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
        ],
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "three_shifted".to_owned(),
        description: "Three Gaussians with shifted means.".to_owned(),
        stress_type: "translation".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.1,
        factors: vec![
          GaussianParams {
            mu: -1.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: 1.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
        ],
        input_grid_domain: (-6.0, 6.0),
        input_grid_n_points: 241,
      },
      GaussianChainTestCase {
        name: "five_identical".to_owned(),
        description: "Five identical Gaussians testing moderate chain.".to_owned(),
        stress_type: "moderate chain".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.15,
        factors: std::iter::repeat_n(
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          5,
        )
        .collect(),
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "ten_identical".to_owned(),
        description: "Ten identical Gaussians testing underflow resistance.".to_owned(),
        stress_type: "underflow risk".to_owned(),
        analytical_caution: "product amplitude decreases rapidly".to_owned(),
        slowness: 0.25,
        factors: std::iter::repeat_n(
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          10,
        )
        .collect(),
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "twenty_identical".to_owned(),
        description: "Twenty identical Gaussians for severe underflow stress.".to_owned(),
        stress_type: "severe underflow".to_owned(),
        analytical_caution: "requires log-scale handling".to_owned(),
        slowness: 0.4,
        factors: std::iter::repeat_n(
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          20,
        )
        .collect(),
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "fifty_identical".to_owned(),
        description: "Fifty identical Gaussians for extreme underflow stress.".to_owned(),
        stress_type: "extreme underflow".to_owned(),
        analytical_caution: "plain multiplication would underflow to zero".to_owned(),
        slowness: 0.7,
        factors: std::iter::repeat_n(
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          50,
        )
        .collect(),
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "five_mixed".to_owned(),
        description: "Five Gaussians with varied parameters.".to_owned(),
        stress_type: "varied params".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.2,
        factors: vec![
          GaussianParams {
            mu: -2.0,
            sigma: 1.5,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: -0.5,
            sigma: 0.8,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: 0.5,
            sigma: 1.2,
            amplitude: 1.0,
          },
          GaussianParams {
            mu: 2.0,
            sigma: 0.9,
            amplitude: 1.0,
          },
        ],
        input_grid_domain: (-8.0, 8.0),
        input_grid_n_points: 321,
      },
      GaussianChainTestCase {
        name: "small_amplitudes_chain".to_owned(),
        description: "Chain of Gaussians with small amplitudes.".to_owned(),
        stress_type: "amplitude underflow".to_owned(),
        analytical_caution: "requires log-scale handling".to_owned(),
        slowness: 0.3,
        factors: std::iter::repeat_n(
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 0.1,
          },
          10,
        )
        .collect(),
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "hundred_identical".to_owned(),
        description: "One hundred identical Gaussians for stress testing.".to_owned(),
        stress_type: "stress test".to_owned(),
        analytical_caution: "extreme precision requirements".to_owned(),
        slowness: 0.9,
        factors: std::iter::repeat_n(
          GaussianParams {
            mu: 0.0,
            sigma: 1.0,
            amplitude: 1.0,
          },
          100,
        )
        .collect(),
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "single_factor".to_owned(),
        description: "Single Gaussian (trivial case, result equals input).".to_owned(),
        stress_type: "edge case".to_owned(),
        analytical_caution: "none".to_owned(),
        slowness: 0.0,
        factors: vec![GaussianParams {
          mu: 0.0,
          sigma: 1.0,
          amplitude: 1.0,
        }],
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
      GaussianChainTestCase {
        name: "empty_factors".to_owned(),
        description: "Empty factors list (edge case, returns identity).".to_owned(),
        stress_type: "edge case".to_owned(),
        analytical_caution: "trivial case".to_owned(),
        slowness: 0.0,
        factors: vec![],
        input_grid_domain: (-5.0, 5.0),
        input_grid_n_points: 201,
      },
    ]
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianChainTestCase {
  pub name: String,
  pub description: String,
  pub stress_type: String,
  pub analytical_caution: String,
  pub slowness: f64,
  pub factors: Vec<GaussianParams>,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl TestCase for GaussianChainTestCase {
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
