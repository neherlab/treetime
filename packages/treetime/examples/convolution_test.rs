use clap::{Parser, ValueEnum};
use eyre::Report;
use serde::{Deserialize, Serialize};
use std::fs;
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use treetime::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use treetime::distribution::reference::convolution_test::framework::framework::ConvolutionTestFramework;
use treetime::distribution::reference::convolution_test::framework::results::TestRunOutcome;
use treetime::distribution::reference::convolution_test::framework::runner::{
  ConvolutionTestRunner, TraitBasedTestRunner,
};
use treetime::distribution::reference::convolution_test::framework::test_case::TestCase;
use treetime::distribution::reference::convolution_test::functions::exponential::ExponentialConvInput;
use treetime::distribution::reference::convolution_test::functions::gaussian::GaussianConvInput;
use treetime::distribution::reference::convolution_test::plots::error_plots::{
  plot_absolute_error, plot_error_histogram, plot_tolerance_metrics,
};
use treetime::distribution::reference::convolution_test::plots::functions_and_convolution::plot_functions_and_convolution;
use treetime::distribution::reference::convolution_test::plots::pointwise_plots::{
  plot_derivative_errors, plot_pointwise_error_profiles,
};
use treetime::distribution::reference::convolution_test::plots::spatial_plots::plot_spatial_profiles;
use treetime::distribution::reference::convolution_test::traits::ConvInput;

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

  for function_type in &args.functions {
    match function_type {
      FunctionType::Gaussian => run_convolution_tests::<GaussianConvInput>(&args)?,
      FunctionType::Exponential => run_convolution_tests::<ExponentialConvInput>(&args)?,
    }
  }

  Ok(())
}

fn run_convolution_tests<I>(args: &Args) -> Result<(), Report>
where
  I: ConvInput,
{
  let function_type_name = I::function_type();

  if args.list_cases {
    I::list_test_cases();
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, function_type_name);
  let runner = TraitBasedTestRunner::<I>::new(Some(args.test_cases.as_str()))?;

  let mut framework = ConvolutionTestFramework::new(runner, output_dir.clone());
  framework.set_algorithms(args.algorithms.clone());

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: {function_type_name}");
    println!("  Algorithms: {:?}", framework.algorithms);
    println!("  Test cases: {} selected", framework.runner.test_cases().len());
    println!("  Output directory: {output_dir}\n");
  }

  let outcomes = framework.run_all_tests()?;

  generate_plot_outputs(&output_dir, &outcomes)?;

  // Generate summary
  let summary = framework.generate_summary(&outcomes);

  // Print summary to console
  framework.print_summary(&summary, &outcomes);

  // Save results
  framework.save_results_json(&outcomes, &summary)?;

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
    }
  }
  Ok(())
}
