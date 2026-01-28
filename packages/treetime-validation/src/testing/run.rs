use crate::algorithms::ConvolutionAlgorithm;
use crate::algorithms::MultiplicationAlgorithm;
use crate::testing::framework::test_case::TestCase;
use crate::testing::runners::convolution::ConvolutionRunner;
use crate::testing::runners::multiplication::ChainMultiplicationRunner;
use crate::testing::runners::multiplication::MultiplicationRunner;
use crate::testing::runners::runner::run_tests_generic;
use crate::testing::test_suites::test_suites::ChainMultiplicationTestSuite;
use crate::testing::test_suites::test_suites::ConvolutionTestSuite;
use crate::testing::test_suites::test_suites::MultiplicationTestSuite;
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
  #[arg(long, value_delimiter = ',', default_values_t = TestSuiteName::all(), help_heading = "Test Selection")]
  pub test_suites: Vec<TestSuiteName>,

  #[arg(long, default_value = "all", help_heading = "Test Selection")]
  pub test_cases: String,

  #[arg(
    long,
    default_value_t = 0.5,
    help_heading = "Test Selection",
    help = "Run test cases with slowness <= this threshold (0.0-1.0)"
  )]
  pub slowness: f64,

  #[arg(long, value_delimiter = ',', default_values_t = ConvolutionAlgorithm::all(), help_heading = "Algorithms")]
  pub conv_algorithms: Vec<ConvolutionAlgorithm>,

  #[arg(long, value_delimiter = ',', default_values_t = MultiplicationAlgorithm::all(), help_heading = "Algorithms")]
  pub mult_algorithms: Vec<MultiplicationAlgorithm>,

  #[arg(long, default_value = "tmp/testing", help_heading = "Output")]
  pub output_dir: String,

  #[arg(
    long,
    help_heading = "Output",
    help = "Show detailed output including per-test metrics"
  )]
  pub verbose: bool,

  #[arg(
    long,
    help_heading = "Output",
    help = "List available test cases without running tests"
  )]
  pub list_cases: bool,
}

pub fn run_validation_tests() -> Result<(), Report> {
  let mut args = Args::parse();
  args.test_suites = TestSuiteName::expand(&args.test_suites);
  args.conv_algorithms = ConvolutionAlgorithm::expand(&args.conv_algorithms);
  args.mult_algorithms = MultiplicationAlgorithm::expand(&args.mult_algorithms);
  for suite_name in &args.test_suites {
    suite_name.run_tests(&args)?;
  }
  Ok(())
}

pub fn run_convolution_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: ConvolutionTestSuite + Default,
{
  run_tests_generic::<ConvolutionRunner<S>>(args, S::default(), |suite| {
    list_test_cases_generic(&suite.create_test_cases(), suite.test_suite_name());
  })
}

pub fn run_multiplication_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: MultiplicationTestSuite + Default,
{
  run_tests_generic::<MultiplicationRunner<S>>(args, S::default(), |suite| {
    list_test_cases_generic(&suite.create_test_cases(), suite.test_suite_name());
  })
}

pub fn run_chain_multiplication_tests_impl<S>(args: &Args) -> Result<(), Report>
where
  S: ChainMultiplicationTestSuite + Default,
{
  run_tests_generic::<ChainMultiplicationRunner<S>>(args, S::default(), |suite| {
    list_test_cases_generic(&suite.create_test_cases(), suite.test_suite_name());
  })
}

fn list_test_cases_generic<T: TestCase>(cases: &[T], suite_name: &str) {
  println!("Available {suite_name} test cases:");
  for case in cases {
    println!(
      "  - {} (slowness: {}) : {}",
      case.name(),
      case.slowness(),
      case.description()
    );
  }
}
