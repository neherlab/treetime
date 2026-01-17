use crate::testing::framework::test_case::TestCase;
use crate::testing::run::{
  Args, run_chain_multiplication_tests_impl, run_convolution_tests_impl, run_multiplication_tests_impl,
};
use crate::testing::test_suites::exponential::ExponentialTestSuite;
use crate::testing::test_suites::gaussian::GaussianTestSuite;
use crate::testing::test_suites::gaussian_chain_multiplication::GaussianChainMultiplicationTestSuite;
use crate::testing::test_suites::gaussian_exponential::GaussianExponentialTestSuite;
use crate::testing::test_suites::gaussian_multiplication::GaussianMultiplicationTestSuite;
use clap::ValueEnum;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use treetime_ops::ScaledArray;
use treetime_utils::make_error;

#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Display, ValueEnum, Serialize, Deserialize, EnumIter)]
#[serde(rename_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum TestSuiteName {
  All,
  #[default]
  Gaussian,
  Exponential,
  GaussianExponential,
  GaussianMultiplication,
  GaussianChainMultiplication,
}

impl TestSuiteName {
  pub fn all() -> Vec<Self> {
    Self::iter().filter(|s| *s != Self::All).collect()
  }

  pub fn convolution_suites() -> Vec<Self> {
    vec![Self::Gaussian, Self::Exponential, Self::GaussianExponential]
  }

  pub fn multiplication_suites() -> Vec<Self> {
    vec![Self::GaussianMultiplication, Self::GaussianChainMultiplication]
  }

  pub fn expand(suites: &[Self]) -> Vec<Self> {
    if suites.contains(&Self::All) {
      Self::all()
    } else {
      suites.to_vec()
    }
  }

  pub fn is_convolution(&self) -> bool {
    matches!(self, Self::Gaussian | Self::Exponential | Self::GaussianExponential)
  }

  pub fn is_multiplication(&self) -> bool {
    matches!(self, Self::GaussianMultiplication | Self::GaussianChainMultiplication)
  }

  /// Run tests for this suite.
  ///
  /// Returns error if called on `All` meta-variant (use `expand()` first).
  pub fn run_tests(&self, args: &Args) -> Result<(), Report> {
    match self {
      Self::All => make_error!("Cannot run All meta-variant; use expand() first"),
      Self::Gaussian => run_convolution_tests_impl::<GaussianTestSuite>(args),
      Self::Exponential => run_convolution_tests_impl::<ExponentialTestSuite>(args),
      Self::GaussianExponential => run_convolution_tests_impl::<GaussianExponentialTestSuite>(args),
      Self::GaussianMultiplication => run_multiplication_tests_impl::<GaussianMultiplicationTestSuite>(args),
      Self::GaussianChainMultiplication => {
        run_chain_multiplication_tests_impl::<GaussianChainMultiplicationTestSuite>(args)
      },
    }
  }
}

pub trait TestSuite: Send + Sync {
  type TestCase: TestCase;

  fn test_suite_name(&self) -> &'static str;

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn analytical_convolution(&self, test_case: &Self::TestCase, eval_grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn create_test_cases(&self) -> Vec<Self::TestCase>;
}

pub trait MultiplicationTestSuite: Send + Sync {
  type TestCase: TestCase;

  fn test_suite_name(&self) -> &'static str;

  fn create_f(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn create_g(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Array1<f64>, Report>;

  fn analytical_multiplication(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<ScaledArray, Report>;

  fn create_test_cases(&self) -> Vec<Self::TestCase>;
}

pub trait ChainMultiplicationTestSuite: Send + Sync {
  type TestCase: TestCase;

  fn test_suite_name(&self) -> &'static str;

  fn create_factors(&self, test_case: &Self::TestCase, grid: &Array1<f64>) -> Result<Vec<Array1<f64>>, Report>;

  fn analytical_chain_multiplication(
    &self,
    test_case: &Self::TestCase,
    grid: &Array1<f64>,
  ) -> Result<ScaledArray, Report>;

  fn create_test_cases(&self) -> Vec<Self::TestCase>;
}
