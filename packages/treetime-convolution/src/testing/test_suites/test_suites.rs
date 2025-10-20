use crate::testing::framework::test_case::TestCase;
use crate::testing::run::{Args, run_convolution_tests_impl};
use crate::testing::test_suites::exponential::ExponentialTestSuite;
use crate::testing::test_suites::gaussian::GaussianTestSuite;
use crate::testing::test_suites::gaussian_exponential::GaussianExponentialTestSuite;
use clap::ValueEnum;
use eyre::Report;
use ndarray::Array1;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};

#[derive(Copy, Clone, Debug, Default, Display, ValueEnum, Serialize, Deserialize, EnumIter)]
#[serde(rename_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum TestSuiteName {
  #[default]
  Gaussian,
  Exponential,
  GaussianExponential,
}

impl TestSuiteName {
  pub fn all() -> Vec<Self> {
    Self::iter().collect()
  }

  pub fn run_tests(&self, args: &Args) -> Result<(), Report> {
    match self {
      Self::Gaussian => run_convolution_tests_impl::<GaussianTestSuite>(args),
      Self::Exponential => run_convolution_tests_impl::<ExponentialTestSuite>(args),
      Self::GaussianExponential => run_convolution_tests_impl::<GaussianExponentialTestSuite>(args),
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
