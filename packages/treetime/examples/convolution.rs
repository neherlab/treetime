use clap::Parser;
use ndarray::Array1;
use plotters::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use treetime::distribution::reference::convolution::{ndarray_convolve, riemann_convolve};
use treetime::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use treetime::distribution::reference::convolution_test::gaussian::analytical::{gaussian_f, gaussian_g};
use treetime::io::json::{JsonPretty, json_write_file, json_write_str};

#[derive(Parser, Clone, Serialize, Deserialize)]
#[command(
  name = "convolution-functions-test",
  about = "Enhanced convolution functions test with metrics and plots",
  version
)]
struct Args {
  /// Algorithm to test
  #[arg(long, default_value = "riemann")]
  algorithm: String,

  /// Standard deviation of first Gaussian f(x)
  #[arg(long, default_value = "1.0")]
  sigma_f: f64,

  /// Standard deviation of second Gaussian g(x)
  #[arg(long, default_value = "1.0")]
  sigma_g: f64,

  /// Mean offset of second Gaussian g(x)
  #[arg(long, default_value = "0.0")]
  mu: f64,

  /// Minimum x value for grid
  #[arg(long, default_value = "-3.0")]
  x_min: f64,

  /// Maximum x value for grid
  #[arg(long, default_value = "3.0")]
  x_max: f64,

  /// Grid spacing
  #[arg(long, default_value = "0.2")]
  dx: f64,

  /// Output directory for plots and results
  #[arg(long, default_value = "tmp/conv_test")]
  output_dir: String,

  /// Range around peak for detailed comparison (±n points)
  #[arg(long, default_value = "5")]
  peak_comparison_range: usize,
}

#[derive(Serialize, Deserialize)]
struct ConvolutionResults {
  parameters: Args,
  function_properties: FunctionProperties,
  grid_info: GridInfo,
  domain_agreement_metrics: DomainAgreementMetrics,
  peak_comparison: Vec<PeakComparisonPoint>,
}

#[derive(Serialize, Deserialize)]
struct FunctionProperties {
  f_variance: f64,
  g_mean: f64,
  g_variance: f64,
  expected_peak_x: f64,
  expected_variance: f64,
  expected_peak_value: f64,
}

#[derive(Serialize, Deserialize)]
struct GridInfo {
  x_min: f64,
  x_max: f64,
  n_points: usize,
  dx: f64,
  algorithm: String,
}

#[derive(Serialize, Deserialize)]
struct PeakComparisonPoint {
  x: f64,
  actual: f64,
  expected: f64,
  absolute_error: f64,
  relative_error_percent: f64,
}

fn main() -> eyre::Result<()> {
  let args = Args::parse();

  println!("=== Enhanced Convolution Functions Test ===\n");

  println!("Args:\n{}", json_write_str(&args, JsonPretty(true))?);

  // Create evaluation grid
  let n_points = ((args.x_max - args.x_min) / args.dx + 1.0).round() as usize;
  let x_grid = Array1::from_iter((0..n_points).map(|i| args.x_min + i as f64 * args.dx));

  // Calculate expected properties
  let expected_peak_x = args.mu;
  let expected_sigma_combined = args.sigma_f.hypot(args.sigma_g);
  let expected_peak_value = 1.0 / (expected_sigma_combined * (2.0 * PI).sqrt());

  let function_properties = FunctionProperties {
    f_variance: args.sigma_f.powi(2),
    g_mean: args.mu,
    g_variance: args.sigma_g.powi(2),
    expected_peak_x,
    expected_variance: expected_sigma_combined.powi(2),
    expected_peak_value,
  };

  println!(
    "Function properties:\n{}",
    json_write_str(&function_properties, JsonPretty(true))?
  );

  let grid_info = GridInfo {
    x_min: args.x_min,
    x_max: args.x_max,
    n_points,
    dx: args.dx,
    algorithm: args.algorithm.clone(),
  };

  println!("Grid info:\n{}", json_write_str(&grid_info, JsonPretty(true))?);

  // Create test functions using grid
  let f = gaussian_f(args.sigma_f, (args.x_min, args.x_max), args.dx)?;
  let g = gaussian_g(args.sigma_g, args.mu, (args.x_min, args.x_max), args.dx)?;

  println!("Testing {} algorithm", args.algorithm);

  let result = match args.algorithm.as_str() {
    "riemann" => riemann_convolve(&f, &g, &x_grid)?,
    "ndarray" => ndarray_convolve(&f, &g, &x_grid)?,
    _ => {
      eprintln!("Unknown algorithm: {}", args.algorithm);
      return Ok(());
    },
  };

  // Compute analytical expected result
  let expected: Array1<f64> = x_grid.mapv(|x| {
    let var_sum = args.sigma_g.powi(2) + args.sigma_f.powi(2);
    (-0.5_f64 * (x - args.mu).powi(2) / var_sum).exp() / (var_sum * 2.0_f64 * PI).sqrt()
  });

  // Comprehensive domain-wide agreement metrics
  let domain_metrics = compute_domain_agreement_metrics(&x_grid, result.y(), &expected)?;

  // Generate plots
  std::fs::create_dir_all(&args.output_dir)?;
  plot_input_functions(&f, &g, &args)?;
  plot_convolution_results(&x_grid, result.y(), &expected, &args)?;
  plot_error_analysis(&x_grid, result.y(), &expected, &args)?;

  // Compare values around peak
  println!("\nComparison around peak (μ = {}):", args.mu);
  println!(
    "{:>8} {:>12} {:>12} {:>12} {:>10}",
    "x", "Actual", "Expected", "Difference", "Rel Error"
  );

  // Find peak region on x_grid
  let peak_x = args.mu;
  let mut peak_idx = 0;
  for (i, &x) in x_grid.iter().enumerate() {
    if (x - peak_x).abs() < (x_grid[peak_idx] - peak_x).abs() {
      peak_idx = i;
    }
  }

  println!("{}", "-".repeat(66));

  let mut peak_comparison_points = Vec::new();

  for i in
    (peak_idx.saturating_sub(args.peak_comparison_range))..(peak_idx + args.peak_comparison_range + 1).min(x_grid.len())
  {
    let x_val = x_grid[i];
    let actual_val = result.y()[i];
    let expected_val = expected[i];
    let diff = (actual_val - expected_val).abs();
    let rel_err = if expected_val != 0.0 {
      diff / expected_val * 100.0
    } else {
      0.0
    };

    peak_comparison_points.push(PeakComparisonPoint {
      x: x_val,
      actual: actual_val,
      expected: expected_val,
      absolute_error: diff,
      relative_error_percent: rel_err,
    });

    let diff_str = if diff < 1e-12 {
      if diff == 0.0 {
        "0".to_owned()
      } else {
        format!("{diff:.2e}")
      }
    } else {
      format!("{diff:.3e}")
    };

    let rel_err_str = if rel_err < 1e-6 {
      if rel_err == 0.0 {
        "0".to_owned()
      } else {
        format!("{rel_err:.2e}")
      }
    } else {
      format!("{rel_err:.2}")
    };

    println!("{x_val:>8.2} {actual_val:>12.6} {expected_val:>12.6} {diff_str:>12} {rel_err_str:>9}%");
  }

  // Create and save results structure
  let results = ConvolutionResults {
    parameters: args.clone(),
    function_properties,
    grid_info,
    domain_agreement_metrics: domain_metrics,
    peak_comparison: peak_comparison_points,
  };

  let results_path = format!("{}/results.json", args.output_dir);
  json_write_file(&results_path, &results, JsonPretty(true))?;
  println!("\nSaved results to {results_path}");

  println!("\nMax value: {:.6}", result.max_value());
  println!("Test completed successfully!");

  Ok(())
}

fn compute_domain_agreement_metrics(
  x_grid: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
) -> eyre::Result<DomainAgreementMetrics> {
  let metrics = DomainAgreementMetrics::new(x_grid, actual, expected)?;
  println!("\n{metrics}");
  Ok(metrics)
}

fn plot_input_functions(
  f: &treetime::distribution::reference::grid_fn::GridFn,
  g: &treetime::distribution::reference::grid_fn::GridFn,
  args: &Args,
) -> eyre::Result<()> {
  let output_path = format!("{}/input_functions.svg", args.output_dir);
  let root = SVGBackend::new(&output_path, (800, 600)).into_drawing_area();
  root.fill(&WHITE)?;

  let max_val = f.y().iter().chain(g.y().iter()).fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption("Input Gaussian Functions", ("Arial", 24))
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(50)
    .build_cartesian_2d(args.x_min..args.x_max, 0.0..max_val * 1.1)?;

  chart.configure_mesh().draw()?;

  let f_data: Vec<(f64, f64)> = f.x().iter().zip(f.y().iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(f_data, BLUE.stroke_width(2)))?
    .label(format!("f(x): N(0, {}²)", args.sigma_f))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], BLUE));

  let g_data: Vec<(f64, f64)> = g.x().iter().zip(g.y().iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(g_data, RED.stroke_width(2)))?
    .label(format!("g(x): N({}, {}²)", args.mu, args.sigma_g))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], RED));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved input functions plot to {output_path}");

  Ok(())
}

fn plot_convolution_results(
  x_grid: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  args: &Args,
) -> eyre::Result<()> {
  let output_path = format!("{}/convolution_results.svg", args.output_dir);
  let root = SVGBackend::new(&output_path, (800, 600)).into_drawing_area();
  root.fill(&WHITE)?;

  let max_val = actual.iter().chain(expected.iter()).fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!("Convolution Results: (f * g)(x) [{}]", args.algorithm),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(50)
    .build_cartesian_2d(args.x_min..args.x_max, 0.0..max_val * 1.1)?;

  chart.configure_mesh().draw()?;

  let actual_data: Vec<(f64, f64)> = x_grid.iter().zip(actual.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(actual_data, BLUE.stroke_width(2)))?
    .label(format!("Actual ({})", args.algorithm))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], BLUE));

  let expected_data: Vec<(f64, f64)> = x_grid.iter().zip(expected.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(expected_data, RED.stroke_width(3)))?
    .label("Expected")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], RED));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved convolution results plot to {output_path}");

  Ok(())
}

fn plot_error_analysis(
  x_grid: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  args: &Args,
) -> eyre::Result<()> {
  let output_path = format!("{}/error_analysis.svg", args.output_dir);
  let root = SVGBackend::new(&output_path, (800, 600)).into_drawing_area();
  root.fill(&WHITE)?;

  let absolute_errors: Vec<f64> = actual
    .iter()
    .zip(expected.iter())
    .map(|(&act, &exp)| (act - exp).abs())
    .collect();

  let max_error = absolute_errors.iter().fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption(
      format!("Absolute Error: |Actual - Expected| [{}]", args.algorithm),
      ("Arial", 24),
    )
    .margin(20)
    .x_label_area_size(40)
    .y_label_area_size(50)
    .build_cartesian_2d(args.x_min..args.x_max, 0.0..max_error * 1.1)?;

  chart.configure_mesh().draw()?;

  let error_data: Vec<(f64, f64)> = x_grid
    .iter()
    .zip(absolute_errors.iter())
    .map(|(&x, y)| (x, *y))
    .collect();
  chart
    .draw_series(LineSeries::new(error_data, GREEN.stroke_width(2)))?
    .label("Absolute Error")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 10, y)], GREEN));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved error analysis plot to {output_path}");

  Ok(())
}
