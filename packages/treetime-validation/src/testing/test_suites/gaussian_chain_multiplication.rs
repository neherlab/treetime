#![allow(clippy::many_single_char_names)]
use crate::testing::framework::test_case::{TestCase, TestCaseBase};
use crate::testing::test_suites::test_suites::ChainMultiplicationTestSuite;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use treetime_analytical::validation::cases::gaussian_chain_multiplication::{
  GaussianChainMultiplicationTestCase as AnalyticalCase, get_gaussian_chain_multiplication_cases,
};
use treetime_analytical::{GaussianParams, gaussian_product};
use treetime_ops::ScaledArray;

#[derive(Default)]
pub struct GaussianChainMultiplicationTestSuite;

impl ChainMultiplicationTestSuite for GaussianChainMultiplicationTestSuite {
  type TestCase = GaussianChainTestCase;

  fn test_suite_name(&self) -> &'static str {
    "mult-gaussian-chain"
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
  ) -> Result<ScaledArray, Report> {
    Ok(gaussian_product(&test_case.factors, grid))
  }

  fn create_test_cases(&self) -> Vec<Self::TestCase> {
    get_gaussian_chain_multiplication_cases()
      .into_iter()
      .map(GaussianChainTestCase::from)
      .collect()
  }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaussianChainTestCase {
  #[serde(flatten)]
  pub base: TestCaseBase,
  pub factors: Vec<GaussianParams>,
  pub input_grid_domain: (f64, f64),
  pub input_grid_n_points: usize,
}

impl From<AnalyticalCase> for GaussianChainTestCase {
  fn from(case: AnalyticalCase) -> Self {
    Self {
      base: TestCaseBase {
        name: case.name.to_owned(),
        description: case.description.to_owned(),
        stress_type: case.stress_type.to_owned(),
        analytical_caution: case.analytical_caution.to_owned(),
        slowness: case.slowness,
      },
      factors: case.factors,
      input_grid_domain: case.input_grid_domain,
      input_grid_n_points: case.input_grid_n_points,
    }
  }
}

impl TestCase for GaussianChainTestCase {
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
