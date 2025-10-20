use clap::{Parser, ValueEnum};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fs;
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use treetime::distribution::reference::convolution_test::framework::{
  ConvolutionTestRunner, TestCase, TestResult, TestRunOutcome, TraitBasedTestRunner,
};
use treetime::distribution::reference::convolution_test::functions::exponential::conv_input::ExponentialConvInput;
use treetime::distribution::reference::convolution_test::functions::exponential::flat_result::ExponentialFlatResult;
use treetime::distribution::reference::convolution_test::functions::gaussian::conv_input::GaussianConvInput;
use treetime::distribution::reference::convolution_test::functions::gaussian::flat_result::GaussianFlatResult;
use treetime::distribution::reference::convolution_test::plots::{
  plot_absolute_error, plot_derivative_errors, plot_error_histogram, plot_functions_and_convolution,
  plot_pointwise_error_profiles, plot_spatial_profiles, plot_tolerance_metrics,
};
use treetime::distribution::reference::convolution_test::traits::ConvInput;
use treetime::distribution::reference::convolution_test::{
  algorithms::ConvolutionAlgorithm,
  framework::ConvolutionTestFramework,
  output::{TestOutputWriter, ToFlatResult},
};

#[derive(Parser, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[command(
  name = "convolution-test",
  about = "Comprehensive convolution accuracy and performance test framework for all function types",
  version
)]
struct Args {
  /// Function types to test
  #[arg(long, default_values_t = FunctionType::all())]
  functions: Vec<FunctionType>,

  /// Output directory for results files
  #[arg(long, default_value = "tmp/convolution_test")]
  output_dir: String,

  /// Algorithms to test
  #[arg(long, default_values_t = ConvolutionAlgorithm::all())]
  algorithms: Vec<ConvolutionAlgorithm>,

  /// Run only specific test cases (comma-separated names, or "all" for all)
  #[arg(long, default_value = "all")]
  test_cases: String,

  /// Enable detailed progress output
  #[arg(long)]
  verbose: bool,

  /// List available test cases for the specified function types
  #[arg(long)]
  list_cases: bool,
}

#[derive(Copy, Clone, Debug, Default, Display, ValueEnum, Serialize, Deserialize, EnumIter)]
#[serde(rename_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
enum FunctionType {
  #[default]
  Gaussian,
  Exponential,
}

impl FunctionType {
  /// Get all available function types
  pub fn all() -> Vec<Self> {
    Self::iter().collect()
  }
}

fn main() -> Result<(), Report> {
  let args = Args::parse();

  if args.list_cases {
    for function_type in &args.functions {
      match function_type {
        FunctionType::Gaussian => GaussianConvInput::new().list_test_cases(),
        FunctionType::Exponential => ExponentialConvInput::new().list_test_cases(),
      }
    }
    return Ok(());
  }

  for function_type in &args.functions {
    match function_type {
      FunctionType::Gaussian => {
        dispatch_function_test::<GaussianConvInput, GaussianFlatResult>(&args)?;
      },
      FunctionType::Exponential => {
        dispatch_function_test::<ExponentialConvInput, ExponentialFlatResult>(&args)?;
      },
    }
  }

  Ok(())
}

fn dispatch_function_test<I, F>(args: &Args) -> Result<(), Report>
where
  I: ConvInput,
  TestResult<I::TestCase>: ToFlatResult<FlatResult = F>,
  F: Serialize,
{
  let input = match args.test_cases.as_str() {
    "all" => I::new(),
    test_cases_str => {
      let filtered = I::new().filter_test_cases(test_cases_str)?;
      I::with_test_cases(filtered)
    },
  };

  let function_type_name = input.function_type();
  let function_output_dir = format!("{}/{}", args.output_dir, function_type_name);
  let runner = TraitBasedTestRunner::new(input);
  run_convolution_tests(args, runner, function_type_name, &function_output_dir)
}

fn run_convolution_tests<R, T, F>(
  args: &Args,
  runner: R,
  function_type_name: &str,
  output_dir: &str,
) -> Result<(), Report>
where
  R: ConvolutionTestRunner<T>,
  T: TestCase,
  TestResult<T>: ToFlatResult<FlatResult = F>,
  F: Serialize,
{
  let mut framework = ConvolutionTestFramework::new(runner, output_dir.to_owned());
  framework.set_algorithms(args.algorithms.clone());

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: {function_type_name}");
    println!("  Algorithms: {:?}", framework.algorithms);
    println!("  Test cases: {} selected", framework.runner.test_cases().len());
    println!("  Output directory: {output_dir}\n");
  }

  let outcomes = framework.run_all_tests()?;

  generate_plot_outputs(output_dir, &outcomes)?;

  // Generate summary
  let summary = framework.generate_summary(&outcomes);

  // Print summary to console
  framework.print_summary(&summary, &outcomes);

  // Save results
  framework.save_results_json(&outcomes, &summary)?;

  // Save TSV with function-specific columns
  let flat_results: Vec<_> = outcomes
    .iter()
    .filter_map(|outcome| match outcome {
      TestRunOutcome::Success(result) => Some(result.to_flat_result()),
      TestRunOutcome::Failure(_) => None,
    })
    .collect();
  let output_writer = TestOutputWriter::new(output_dir.to_owned());
  output_writer.save_results_tsv(&flat_results)?;

  println!("{function_type_name} convolution test framework completed successfully!");
  println!("Check {output_dir} for detailed results.");

  Ok(())
}

fn generate_plot_outputs<T>(output_dir: &str, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  for outcome in outcomes {
    if let TestRunOutcome::Success(result) = outcome {
      let algorithm_dir = format!(
        "{}/{}/{}",
        output_dir,
        result.test_case.name(),
        result.algorithm.to_string().to_lowercase()
      );
      fs::create_dir_all(&algorithm_dir)?;

      plot_functions_and_convolution(result, &algorithm_dir)?;
      plot_absolute_error(result, &algorithm_dir)?;
      plot_tolerance_metrics(result, &algorithm_dir)?;
      plot_pointwise_error_profiles(result, &algorithm_dir)?;
      plot_derivative_errors(result, &algorithm_dir)?;
      plot_spatial_profiles(result, &algorithm_dir)?;
      plot_error_histogram(result, &algorithm_dir)?;

      let output_writer = TestOutputWriter::new(output_dir.to_owned());
      output_writer.save_pointwise_arrays(result, &algorithm_dir)?;
      output_writer.save_error_histogram(result, &algorithm_dir)?;
    }
  }
  Ok(())
}
