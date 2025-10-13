use clap::{Parser, ValueEnum};
use eyre::Report;
use itertools::izip;
use ndarray::Array1;
use plotters::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs;
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

  generate_plot_outputs(args, &outcomes)?;

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

fn generate_plot_outputs<T>(args: &Args, outcomes: &[TestRunOutcome<T>]) -> Result<(), Report>
where
  T: TestCase,
{
  for outcome in outcomes {
    if let TestRunOutcome::Success(result) = outcome {
      let algorithm_dir = format!(
        "{}/{}/{}",
        args.output_dir,
        result.test_case.name(),
        result.algorithm.to_string().to_lowercase()
      );
      fs::create_dir_all(&algorithm_dir)?;
      plot_functions_and_convolution(result, &algorithm_dir)?;
      plot_absolute_error(result, &algorithm_dir)?;
      plot_tolerance_metrics(result, &algorithm_dir)?;
    }
  }
  Ok(())
}

fn plot_functions_and_convolution<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let output_path = format!("{output_dir}/functions_convolution.svg");
  let root = SVGBackend::new(&output_path, (960, 640)).into_drawing_area();
  root.fill(&WHITE)?;
  let (x_raw_min, x_raw_max) = combined_range(&[&result.f_x_values, &result.g_x_values, &result.evaluation_grid]);
  let (x_min, x_max) = expand_range(x_raw_min, x_raw_max);
  let (y_raw_min, y_raw_max) = combined_range(&[
    &result.f_y_values,
    &result.g_y_values,
    &result.actual_values,
    &result.expected_values,
  ]);
  let (y_min, y_max) = expand_range(y_raw_min, y_raw_max);
  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!(
        "{} | {} | Input and Convolution",
        result.test_case.name(),
        result.algorithm
      ),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(60)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
  chart.configure_mesh().draw()?;
  chart
    .draw_series(LineSeries::new(
      result
        .f_x_values
        .iter()
        .zip(result.f_y_values.iter())
        .map(|(&x, &y)| (x, y)),
      BLUE.stroke_width(2),
    ))?
    .label("f(x)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], BLUE.stroke_width(2)));
  chart
    .draw_series(LineSeries::new(
      result
        .g_x_values
        .iter()
        .zip(result.g_y_values.iter())
        .map(|(&x, &y)| (x, y)),
      GREEN.stroke_width(2),
    ))?
    .label("g(x)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], GREEN.stroke_width(2)));
  chart
    .draw_series(LineSeries::new(
      result
        .evaluation_grid
        .iter()
        .zip(result.actual_values.iter())
        .map(|(&x, &y)| (x, y)),
      MAGENTA.stroke_width(3),
    ))?
    .label(format!("Actual ({})", result.algorithm))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], MAGENTA.stroke_width(3)));
  chart
    .draw_series(LineSeries::new(
      result
        .evaluation_grid
        .iter()
        .zip(result.expected_values.iter())
        .map(|(&x, &y)| (x, y)),
      RED.stroke_width(3),
    ))?
    .label("Expected")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], RED.stroke_width(3)));
  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE.mix(0.8))
    .draw()?;
  root.present()?;
  Ok(())
}

fn plot_absolute_error<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let mut error_points = Vec::with_capacity(result.evaluation_grid.len());
  let mut max_error = 0.0_f64;
  for (&x, &actual, &expected) in izip!(
    result.evaluation_grid.iter(),
    result.actual_values.iter(),
    result.expected_values.iter()
  ) {
    let value = (actual - expected).abs();
    if value > max_error {
      max_error = value;
    }
    error_points.push((x, value));
  }
  let (x_raw_min, x_raw_max) = combined_range(&[&result.evaluation_grid]);
  let (x_min, x_max) = expand_range(x_raw_min, x_raw_max);
  let (mut y_min, mut y_max) = expand_range(0.0, max_error);
  if y_min < 0.0 {
    y_min = 0.0;
  }
  if y_max <= y_min {
    y_max = y_min + 1.0;
  }
  let output_path = format!("{output_dir}/absolute_error.svg");
  let root = SVGBackend::new(&output_path, (960, 480)).into_drawing_area();
  root.fill(&WHITE)?;
  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!("{} | {} | Absolute Error", result.test_case.name(), result.algorithm),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(60)
    .build_cartesian_2d(x_min..x_max, y_min..y_max)?;
  chart.configure_mesh().y_desc("|actual - expected|").draw()?;
  chart
    .draw_series(LineSeries::new(error_points.iter().copied(), RED.stroke_width(2)))?
    .label("Absolute error")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 12, y)], RED.stroke_width(2)));
  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE.mix(0.8))
    .draw()?;
  root.present()?;
  Ok(())
}

fn plot_tolerance_metrics<T>(result: &TestResult<T>, output_dir: &str) -> Result<(), Report>
where
  T: TestCase,
{
  let abs = result.metrics.abs_tolerance_fractions;
  let rel = result.metrics.rel_tolerance_fractions;
  let output_path = format!("{output_dir}/tolerance_metrics.svg");
  let root = SVGBackend::new(&output_path, (800, 600)).into_drawing_area();
  root.fill(&WHITE)?;
  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!(
        "{} | {} | Tolerance Fractions",
        result.test_case.name(),
        result.algorithm
      ),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(50)
    .y_label_area_size(60)
    .build_cartesian_2d(-0.2..3.0, 0.0..1.05)?;
  chart
    .configure_mesh()
    .x_labels(3)
    .y_labels(6)
    .x_desc("Tolerance level")
    .y_desc("Fraction within threshold")
    .x_label_formatter(&|value| tolerance_label(*value))
    .y_label_formatter(&|value| format!("{value:.2}"))
    .draw()?;
  chart
    .draw_series((0..3).map(|i| {
      let x0 = i as f64 + 0.05;
      let x1 = x0 + 0.35;
      Rectangle::new([(x0, 0.0), (x1, abs[i])], BLUE.mix(0.6).filled())
    }))?
    .label("Absolute")
    .legend(|(x, y)| Rectangle::new([(x - 6, y - 6), (x + 6, y + 6)], BLUE.mix(0.6).filled()));
  chart
    .draw_series((0..3).map(|i| {
      let x0 = i as f64 + 0.5;
      let x1 = x0 + 0.35;
      Rectangle::new([(x0, 0.0), (x1, rel[i])], RED.mix(0.6).filled())
    }))?
    .label("Relative")
    .legend(|(x, y)| Rectangle::new([(x - 6, y - 6), (x + 6, y + 6)], RED.mix(0.6).filled()));
  chart
    .configure_series_labels()
    .border_style(BLACK)
    .background_style(WHITE.mix(0.8))
    .draw()?;
  root.present()?;
  Ok(())
}

fn combined_range(arrays: &[&Array1<f64>]) -> (f64, f64) {
  let mut min_value = f64::INFINITY;
  let mut max_value = f64::NEG_INFINITY;
  for array in arrays {
    for &value in *array {
      if value < min_value {
        min_value = value;
      }
      if value > max_value {
        max_value = value;
      }
    }
  }
  if min_value.is_finite() && max_value.is_finite() {
    (min_value, max_value)
  } else {
    (0.0, 1.0)
  }
}

fn expand_range(min_value: f64, max_value: f64) -> (f64, f64) {
  if !min_value.is_finite() || !max_value.is_finite() {
    return (-1.0, 1.0);
  }
  if max_value <= min_value {
    let base = if min_value.abs() > 1.0 { min_value.abs() } else { 1.0 };
    return (min_value - base, max_value + base);
  }
  let span = max_value - min_value;
  let padding = if span * 0.05 < 1e-9 { 1e-9 } else { span * 0.05 };
  (min_value - padding, max_value + padding)
}

fn tolerance_label(value: f64) -> String {
  let idx = (value + 0.5).floor() as i32;
  match idx {
    0 => "Strict".to_owned(),
    1 => "Moderate".to_owned(),
    2 => "Loose".to_owned(),
    _ => format!("{value:.1}"),
  }
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
