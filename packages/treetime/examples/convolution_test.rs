use clap::{Parser, ValueEnum};
use eyre::Report;
use serde::{Deserialize, Serialize};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter};
use treetime::distribution::reference::convolution_test::exponential::ExponentialTestCase;
use treetime::distribution::reference::convolution_test::framework::{
  ConvolutionTestRunner, TestCase, TestResult, TestRunOutcome,
};
use treetime::distribution::reference::convolution_test::gaussian::GaussianTestCase;
use treetime::distribution::reference::convolution_test::{
  ConvolutionAlgorithm, ExponentialTestRunner, GaussianTestRunner, GenericConvolutionTestFramework, TestOutputWriter,
  ToFlatResult,
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
        FunctionType::Gaussian => {
          list_test_cases::<GaussianTestRunner, GaussianTestCase>();
        },
        FunctionType::Exponential => {
          list_test_cases::<ExponentialTestRunner, ExponentialTestCase>();
        },
      }
    }
    return Ok(());
  }

  for function_type in &args.functions {
    // Create function-specific output directory
    let function_output_dir = format!("{}/{}", args.output_dir, function_type.to_string().to_lowercase());
    let function_args = Args {
      functions: vec![*function_type],
      output_dir: function_output_dir,
      algorithms: args.algorithms.clone(),
      test_cases: args.test_cases.clone(),
      verbose: args.verbose,
      list_cases: false,
    };

    match function_type {
      FunctionType::Gaussian => {
        run_tests_for_runner::<GaussianTestRunner, GaussianTestCase>(&function_args)?;
      },
      FunctionType::Exponential => {
        run_tests_for_runner::<ExponentialTestRunner, ExponentialTestCase>(&function_args)?;
      },
    }
  }

  Ok(())
}

fn filter_test_cases<R, T>(args: &Args) -> Result<Option<Vec<T>>, Report>
where
  R: TestRunner<T>,
  T: TestCase,
{
  let runner = R::new();
  let all_cases = runner.test_cases();

  if args.test_cases == "all" {
    return Ok(Some(all_cases.to_vec()));
  }

  let requested_cases: Vec<&str> = args.test_cases.split(',').map(|s| s.trim()).collect();
  let filtered_cases = all_cases
    .iter()
    .filter(|case| requested_cases.contains(&case.name()))
    .cloned()
    .collect::<Vec<_>>();

  if filtered_cases.is_empty() {
    eprintln!("No matching test cases found for: {}", args.test_cases);
    eprintln!("Available {} test cases:", R::function_type_name());
    for case in all_cases {
      eprintln!("  - {}", case.name());
    }
    return Ok(None);
  }

  Ok(Some(filtered_cases))
}

fn run_convolution_tests<R, T>(args: &Args, runner: R, function_type_name: &str) -> Result<(), Report>
where
  R: ConvolutionTestRunner<T>,
  T: TestCase,
  TestResult<T>: ToFlatResult,
  <TestResult<T> as ToFlatResult>::FlatResult: Serialize,
{
  // Create framework
  let mut framework = GenericConvolutionTestFramework::new(runner, args.output_dir.clone());
  framework.set_algorithms(args.algorithms.clone());

  if args.verbose {
    println!("Test Configuration:");
    println!("  Function type: {function_type_name}");
    println!("  Algorithms: {:?}", framework.algorithms);
    println!("  Test cases: {} selected", framework.runner.test_cases().len());
    println!("  Output directory: {}\n", args.output_dir);
  }

  // Run all tests
  let outcomes = framework.run_all_tests()?;

  // Generate summary
  let summary = framework.generate_summary(&outcomes);

  // Print summary to console
  framework.print_summary(&summary);

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
  let output_writer = TestOutputWriter::new(args.output_dir.clone());
  output_writer.save_results_tsv(&flat_results)?;

  println!("{function_type_name} convolution test framework completed successfully!");
  println!("Check {} for detailed results.", args.output_dir);

  Ok(())
}

trait TestRunner<T: TestCase>: ConvolutionTestRunner<T> {
  /// Create new test runner with default test cases
  fn new() -> Self;

  /// Create test runner with custom test cases
  fn with_test_cases(test_cases: Vec<T>) -> Self;

  /// Get the function type name
  fn function_type_name() -> &'static str;
}

impl TestRunner<GaussianTestCase> for GaussianTestRunner {
  fn new() -> Self {
    GaussianTestRunner::new()
  }

  fn with_test_cases(test_cases: Vec<GaussianTestCase>) -> Self {
    GaussianTestRunner::with_test_cases(test_cases)
  }

  fn function_type_name() -> &'static str {
    "Gaussian"
  }
}

impl TestRunner<ExponentialTestCase> for ExponentialTestRunner {
  fn new() -> Self {
    ExponentialTestRunner::new()
  }

  fn with_test_cases(test_cases: Vec<ExponentialTestCase>) -> Self {
    ExponentialTestRunner::with_test_cases(test_cases)
  }

  fn function_type_name() -> &'static str {
    "Exponential"
  }
}
fn run_tests_for_runner<R, T>(args: &Args) -> Result<(), Report>
where
  R: TestRunner<T>,
  T: TestCase,
  TestResult<T>: ToFlatResult,
  <TestResult<T> as ToFlatResult>::FlatResult: Serialize,
{
  let filtered_cases = if let Some(cases) = filter_test_cases::<R, T>(args)? {
    cases
  } else {
    return Ok(());
  };

  let runner = if args.test_cases == "all" {
    R::new()
  } else {
    R::with_test_cases(filtered_cases)
  };

  run_convolution_tests(args, runner, R::function_type_name())
}

fn list_test_cases<R, T>()
where
  R: TestRunner<T>,
  T: TestCase,
{
  println!("Available {} test cases:", R::function_type_name());
  for case in R::new().test_cases() {
    println!("  - {} : {}", case.name(), case.description());
  }
}
