use crate::algorithms::ConvolutionAlgorithm;
use crate::algorithms::MultiplicationAlgorithm;
use crate::testing::framework::test_case::TestCase;
use crate::testing::runners::convolution::ConvolutionRunner;
use crate::testing::runners::multiplication::ChainMultiplicationRunner;
use crate::testing::runners::multiplication::MultiplicationRunner;
use crate::testing::runners::runner::run_tests_generic;
use crate::testing::test_suites::test_suites::ChainMultiplicationTestSuite;
use crate::testing::test_suites::test_suites::MultiplicationTestSuite;
use crate::testing::test_suites::test_suites::TestSuite;
use crate::testing::test_suites::test_suites::TestSuiteName;
use clap::Parser;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Parser, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[command(
  name = "validation-test",
  about = "Validation test framework for convolution and multiplication algorithms",
  version
)]
pub struct Args {
  #[arg(long, value_delimiter = ',', default_values_t = TestSuiteName::all())]
  pub test_suites: Vec<TestSuiteName>,

  #[arg(long, default_value = "tmp/testing")]
  pub output_dir: String,

  #[arg(long, value_delimiter = ',', default_values_t = ConvolutionAlgorithm::all())]
  pub algorithms: Vec<ConvolutionAlgorithm>,

  #[arg(long, value_delimiter = ',', default_values_t = MultiplicationAlgorithm::all())]
  pub mult_algorithms: Vec<MultiplicationAlgorithm>,

  #[arg(long, default_value = "all")]
  pub test_cases: String,

  #[arg(long, default_value_t = 0.5)]
  pub slowness: f64,

  #[arg(long)]
  pub verbose: bool,

  #[arg(long)]
  pub list_cases: bool,
}

pub fn run_validation_tests() -> Result<(), Report> {
  let mut args = Args::parse();
  args.test_suites = TestSuiteName::expand(&args.test_suites);
  args.algorithms = ConvolutionAlgorithm::expand(&args.algorithms);
  args.mult_algorithms = MultiplicationAlgorithm::expand(&args.mult_algorithms);
  for suite_name in &args.test_suites {
    suite_name.run_tests(&args)?;
  }
  Ok(())
}

pub fn run_convolution_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: TestSuite + Default,
{
  run_tests_generic::<ConvolutionRunner<S>>(args, S::default(), |suite| list_test_cases(suite))
}

pub fn run_multiplication_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: MultiplicationTestSuite + Default,
{
  run_tests_generic::<MultiplicationRunner<S>>(args, S::default(), list_multiplication_test_cases)
}

pub fn run_chain_multiplication_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: ChainMultiplicationTestSuite + Default,
{
  run_tests_generic::<ChainMultiplicationRunner<S>>(args, S::default(), list_chain_multiplication_test_cases)
}

fn list_test_cases<S: TestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.test_suite_name());
  for case in suite.create_test_cases() {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}

fn list_multiplication_test_cases<S: MultiplicationTestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.test_suite_name());
  for case in suite.create_test_cases() {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}

fn list_chain_multiplication_test_cases<S: ChainMultiplicationTestSuite>(suite: &S) {
  println!("Available {} test cases:", suite.test_suite_name());
  for case in suite.create_test_cases() {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}
