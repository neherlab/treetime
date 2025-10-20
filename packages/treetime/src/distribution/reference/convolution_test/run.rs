use crate::distribution::reference::convolution_test::algorithms::ConvolutionAlgorithm;
use crate::distribution::reference::convolution_test::framework::framework::ConvolutionTestFramework;
use crate::distribution::reference::convolution_test::framework::runner::{
  ConvolutionTestRunner, TraitBasedTestRunner,
};
use crate::distribution::reference::convolution_test::functions::functions::{FunctionType, get_function};
use crate::distribution::reference::convolution_test::plots::plots::generate_plot_outputs;
use crate::distribution::reference::convolution_test::traits::ConvInput;
use clap::Parser;
use eyre::Report;
use serde::{Deserialize, Serialize};

#[derive(Parser, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[command(
  name = "convolution-test",
  about = "Comprehensive convolution accuracy and performance test framework for all function types",
  version
)]
pub struct Args {
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

pub fn run_convolution_tests() -> Result<(), Report> {
  let args = Args::parse();
  for function_type in &args.functions {
    let input = get_function(function_type, &args.test_cases)?;
    run_convolution_tests_impl(input, &args)?;
  }
  Ok(())
}

fn run_convolution_tests_impl<I>(input: I, args: &Args) -> Result<(), Report>
where
  I: ConvInput,
{
  if args.list_cases {
    input.list_test_cases();
    return Ok(());
  }

  let output_dir = format!("{}/{}", args.output_dir, input.function_type());
  let runner = TraitBasedTestRunner::<I>::new(Some(args.test_cases.as_str()))?;

  let mut framework = ConvolutionTestFramework::new(runner, output_dir.clone());
  framework.set_algorithms(args.algorithms.clone());

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: {}", input.function_type());
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

  println!(
    "{} convolution test framework completed successfully!",
    input.function_type()
  );
  println!("Check {output_dir} for detailed results.");

  Ok(())
}
