use clap::{Parser, ValueEnum};
use eyre::Report;
use serde::{Deserialize, Serialize};
use strum_macros::EnumIter;
use treetime::distribution::reference::convolution_test::framework::ConvolutionTestRunner;
use treetime::distribution::reference::convolution_test::{
  ConvolutionAlgorithm, ExponentialFlatResult, ExponentialTestRunner, GaussianFlatResult, GaussianTestRunner,
  GenericConvolutionTestFramework, TestOutputWriter, ToFlatResult,
};

#[derive(Parser, Clone, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[command(
  name = "convolution-test",
  about = "Comprehensive convolution accuracy and performance test framework for all function types",
  version
)]
struct Args {
  /// Function type to test
  #[arg(long, default_value = "gaussian")]
  function: FunctionType,

  /// Output directory for results files
  #[arg(long, default_value = "tmp/convolution_test")]
  output_dir: String,

  /// Algorithms to test
  #[arg(long, default_values = ["riemann", "ndarray"])]
  algorithms: Vec<ConvolutionAlgorithm>,

  /// Run only specific test cases (comma-separated names, or "all" for all)
  #[arg(long, default_value = "all")]
  test_cases: String,

  /// Enable detailed progress output
  #[arg(long)]
  verbose: bool,

  /// List available test cases for the specified function type
  #[arg(long)]
  list_cases: bool,
}

#[derive(Clone, Debug, ValueEnum, Serialize, Deserialize, EnumIter)]
#[serde(rename_all = "kebab-case")]
#[clap(rename_all = "kebab-case")]
enum FunctionType {
  Gaussian,
  Exponential,
}

fn main() -> Result<(), Report> {
  let args = Args::parse();

  match args.function {
    FunctionType::Gaussian => {
      if args.list_cases {
        list_gaussian_test_cases();
        return Ok(());
      }
      run_gaussian_tests(&args)
    },
    FunctionType::Exponential => {
      if args.list_cases {
        list_exponential_test_cases();
        return Ok(());
      }
      run_exponential_tests(&args)
    },
  }
}

fn run_gaussian_tests(args: &Args) -> Result<(), Report> {
  // Create Gaussian test runner
  let runner = GaussianTestRunner::new();

  // Create framework
  let mut framework = GenericConvolutionTestFramework::new(runner, args.output_dir.clone());
  framework.set_algorithms(args.algorithms.clone()); // Filter test cases if requested
  if args.test_cases != "all" {
    let requested_cases: Vec<&str> = args.test_cases.split(',').map(|s| s.trim()).collect();
    let all_cases = framework.runner.test_cases();
    let filtered_cases = all_cases
      .iter()
      .filter(|case| requested_cases.contains(&case.name.as_str()))
      .cloned()
      .collect::<Vec<_>>();

    if filtered_cases.is_empty() {
      eprintln!("No matching test cases found for: {}", args.test_cases);
      eprintln!("Available Gaussian test cases:");
      for case in all_cases {
        eprintln!("  - {}", case.name);
      }
      return Ok(());
    }

    let filtered_runner = GaussianTestRunner::with_test_cases(filtered_cases);
    framework = GenericConvolutionTestFramework::new(filtered_runner, args.output_dir.clone());
    framework.set_algorithms(args.algorithms.clone());
  }

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: Gaussian");
    println!("  Algorithms: {:?}", framework.algorithms);
    println!("  Test cases: {} selected", framework.runner.test_cases().len());
    println!("  Output directory: {}\n", args.output_dir);
  }

  // Run all tests
  let results = framework.run_all_tests()?;

  // Generate summary
  let summary = framework.generate_summary(&results);

  // Print summary to console
  framework.print_summary(&summary);

  // Save results
  framework.save_results_json(&results, &summary)?;

  // Save TSV with Gaussian-specific columns
  let flat_results: Vec<GaussianFlatResult> = results.iter().map(|r| r.to_flat_result()).collect();
  let output_writer = TestOutputWriter::new(args.output_dir.clone());
  output_writer.save_results_tsv(&results, &flat_results)?;

  println!("Gaussian convolution test framework completed successfully!");
  println!("Check {} for detailed results.", args.output_dir);

  Ok(())
}

fn run_exponential_tests(args: &Args) -> Result<(), Report> {
  // Create Exponential test runner
  let runner = ExponentialTestRunner::new();

  // Create framework
  let mut framework = GenericConvolutionTestFramework::new(runner, args.output_dir.clone());
  framework.set_algorithms(args.algorithms.clone()); // Filter test cases if requested
  if args.test_cases != "all" {
    let requested_cases: Vec<&str> = args.test_cases.split(',').map(|s| s.trim()).collect();
    let all_cases = framework.runner.test_cases();
    let filtered_cases = all_cases
      .iter()
      .filter(|case| requested_cases.contains(&case.name.as_str()))
      .cloned()
      .collect::<Vec<_>>();

    if filtered_cases.is_empty() {
      eprintln!("No matching test cases found for: {}", args.test_cases);
      eprintln!("Available Exponential test cases:");
      for case in all_cases {
        eprintln!("  - {}", case.name);
      }
      return Ok(());
    }

    let filtered_runner = ExponentialTestRunner::with_test_cases(filtered_cases);
    framework = GenericConvolutionTestFramework::new(filtered_runner, args.output_dir.clone());
    framework.set_algorithms(args.algorithms.clone());
  }

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: Exponential");
    println!("  Algorithms: {:?}", framework.algorithms);
    println!("  Test cases: {} selected", framework.runner.test_cases().len());
    println!("  Output directory: {}\n", args.output_dir);
  }

  // Run all tests
  let results = framework.run_all_tests()?;

  // Generate summary
  let summary = framework.generate_summary(&results);

  // Print summary to console
  framework.print_summary(&summary);

  // Save results
  framework.save_results_json(&results, &summary)?;

  // Save TSV with Exponential-specific columns
  let flat_results: Vec<ExponentialFlatResult> = results.iter().map(|r| r.to_flat_result()).collect();
  let output_writer = TestOutputWriter::new(args.output_dir.clone());
  output_writer.save_results_tsv(&results, &flat_results)?;

  println!("Exponential convolution test framework completed successfully!");
  println!("Check {} for detailed results.", args.output_dir);

  Ok(())
}

fn list_gaussian_test_cases() {
  let runner = GaussianTestRunner::new();
  println!("Available Gaussian test cases:");
  for case in runner.test_cases() {
    println!("  - {} : {}", case.name, case.description);
  }
}

fn list_exponential_test_cases() {
  let runner = ExponentialTestRunner::new();
  println!("Available Exponential test cases:");
  for case in runner.test_cases() {
    println!("  - {} : {}", case.name, case.description);
  }
}
