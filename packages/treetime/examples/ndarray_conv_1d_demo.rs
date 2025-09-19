use clap::Parser;
use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use plotters::prelude::*;
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use treetime::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use treetime::io::json::{JsonPretty, json_write_file, json_write_str};
use treetime::utils::float_fmt::float_to_digits;

#[derive(Parser, Clone, Serialize, Deserialize)]
#[command(
  name = "ndarray-conv-1d-demo",
  about = "1D Gaussian convolution demo using ndarray-conv",
  version
)]
struct Args {
  /// Standard deviation of first Gaussian f(x)
  #[arg(long, default_value = "1.0")]
  sigma_f: f64,

  /// Standard deviation of second Gaussian g(x)
  #[arg(long, default_value = "1.0")]
  sigma_g: f64,

  /// Mean offset of second Gaussian g(x)
  #[arg(long, default_value = "1.0")]
  mu: f64,

  /// Minimum x value for grid
  #[arg(long, default_value = "-7.0")]
  x_min: f64,

  /// Maximum x value for grid
  #[arg(long, default_value = "9.0")]
  x_max: f64,

  /// Number of grid points
  #[arg(long, default_value = "321")]
  n_points: usize,

  /// Output directory for plots
  #[arg(long, default_value = "tmp/conv_plots")]
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
  conv_len: usize,
  start_idx: usize,
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

  println!("=== Gaussian Convolution Demo with ndarray-conv ===\n");

  println!("Args:\n{}", json_write_str(&args, JsonPretty(true))?);

  // Use a narrower grid centered around the effective convolution domain
  // This allows ndarray-conv to work directly without complex interpolation
  let x_grid = Array1::linspace(args.x_min, args.x_max, args.n_points);
  let dx = (args.x_max - args.x_min) / (args.n_points - 1) as f64;

  // Calculate expected properties
  let expected_peak_x = args.mu;
  let expected_sigma_combined = (args.sigma_f.powi(2) + args.sigma_g.powi(2)).sqrt();
  let expected_peak_value = 1.0 / (expected_sigma_combined * (2.0 * PI).sqrt());

  // Create function properties structure
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

  // Define Gaussian functions on this grid
  let f_vals: Array1<f64> =
    x_grid.mapv(|x| (-0.5_f64 * (x / args.sigma_f).powi(2)).exp() / (args.sigma_f * (2.0_f64 * PI).sqrt()));

  let g_vals: Array1<f64> =
    x_grid.mapv(|x| (-0.5_f64 * ((x - args.mu) / args.sigma_g).powi(2)).exp() / (args.sigma_g * (2.0_f64 * PI).sqrt()));

  // Perform convolution using ndarray-conv with Full mode
  let raw_conv = f_vals.conv(&g_vals, ConvMode::Full, PaddingMode::Zeros)?;

  // Scale by grid spacing to approximate continuous convolution integral
  let conv_scaled = &raw_conv * dx;

  // The Full convolution result spans a wider domain
  // Create the extended domain for the full convolution result
  let conv_len = raw_conv.len();
  let conv_x_min = 2.0 * args.x_min;
  let conv_x_max = 2.0 * args.x_max;
  let conv_x_grid = Array1::linspace(conv_x_min, conv_x_max, conv_len);

  // Extract the portion that corresponds to our original evaluation domain
  // Find indices in conv_x_grid that correspond to our desired x_grid range
  let mut start_idx = 0;
  for (i, &conv_x) in conv_x_grid.iter().enumerate() {
    if conv_x >= args.x_min {
      start_idx = i;
      break;
    }
  }

  // Create grid info
  let grid_info = GridInfo {
    x_min: args.x_min,
    x_max: args.x_max,
    n_points: args.n_points,
    dx,
    conv_len,
    start_idx,
  };

  println!("Grid info:\n{}", json_write_str(&grid_info, JsonPretty(true))?);

  let actual = conv_scaled
    .slice(ndarray::s![start_idx..start_idx + args.n_points])
    .to_owned();

  // Compute expected result on the same grid
  let expected: Array1<f64> = x_grid.mapv(|x| {
    let var_sum = args.sigma_g.powi(2) + args.sigma_f.powi(2);
    (-0.5_f64 * (x - args.mu).powi(2) / var_sum).exp() / (var_sum * 2.0_f64 * PI).sqrt()
  });

  // Comprehensive domain-wide agreement metrics
  let domain_metrics = compute_domain_agreement_metrics(&x_grid, &actual, &expected)?;

  // Generate plots
  plot_input_functions(
    &x_grid,
    &f_vals,
    &g_vals,
    args.sigma_f,
    args.sigma_g,
    args.mu,
    &args.output_dir,
  )?;
  plot_convolution_results(&x_grid, &actual, &expected, &args.output_dir)?;
  plot_error_analysis(&x_grid, &actual, &expected, &args.output_dir)?;

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
    let actual_val = actual[i];
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
      float_to_digits(diff, Some(3), None)
    };

    let rel_err_str = if rel_err < 1e-6 {
      if rel_err == 0.0 {
        "0".to_owned()
      } else {
        format!("{rel_err:.2e}")
      }
    } else {
      float_to_digits(rel_err, Some(3), Some(2))
    };

    println!(
      "{:>8} {:>12} {:>12} {:>12} {:>9}%",
      float_to_digits(x_val, Some(3), Some(2)),
      float_to_digits(actual_val, Some(6), None),
      float_to_digits(expected_val, Some(6), None),
      diff_str,
      rel_err_str
    );
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
  x_grid: &Array1<f64>,
  f_vals: &Array1<f64>,
  g_vals: &Array1<f64>,
  sigma_f: f64,
  sigma_g: f64,
  mu: f64,
  output_dir: &str,
) -> eyre::Result<()> {
  // Plot configuration
  let plot_width = 800;
  let plot_height = 600;
  let plot_margin = 20;
  let x_label_area = 40;
  let y_label_area = 50;
  let line_width = 2;
  let legend_offset = 10;
  let max_val_margin = 1.1;
  let font_size = 24;

  let output_path = format!("{output_dir}/input_functions.svg");
  let root = SVGBackend::new(&output_path, (plot_width, plot_height)).into_drawing_area();
  root.fill(&WHITE)?;

  let max_val = f_vals.iter().chain(g_vals.iter()).fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption("Input Gaussian Functions", ("Arial", font_size))
    .margin(plot_margin)
    .x_label_area_size(x_label_area)
    .y_label_area_size(y_label_area)
    .build_cartesian_2d(-7.0..9.0, 0.0..max_val * max_val_margin)?;

  chart.configure_mesh().draw()?;

  // Plot f(x)
  let f_data: Vec<(f64, f64)> = x_grid.iter().zip(f_vals.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(f_data, BLUE.stroke_width(line_width)))?
    .label(format!("f(x): N(0, {sigma_f}²)"))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], BLUE));

  // Plot g(x)
  let g_data: Vec<(f64, f64)> = x_grid.iter().zip(g_vals.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(g_data, RED.stroke_width(line_width)))?
    .label(format!("g(x): N({mu}, {sigma_g}²)"))
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], RED));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved input functions plot to {output_path}");

  Ok(())
}

fn plot_convolution_results(
  x_grid: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  output_dir: &str,
) -> eyre::Result<()> {
  // Plot configuration
  let plot_width = 800;
  let plot_height = 600;
  let plot_margin = 20;
  let x_label_area = 40;
  let y_label_area = 50;
  let actual_line_width = 2;
  let expected_line_width = 3;
  let legend_offset = 10;
  let max_val_margin = 1.1;
  let font_size = 24;

  let output_path = format!("{output_dir}/convolution_results.svg");
  let root = SVGBackend::new(&output_path, (plot_width, plot_height)).into_drawing_area();
  root.fill(&WHITE)?;

  let max_val = actual.iter().chain(expected.iter()).fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption("Convolution Results: (f * g)(x)", ("Arial", font_size))
    .margin(plot_margin)
    .x_label_area_size(x_label_area)
    .y_label_area_size(y_label_area)
    .build_cartesian_2d(-7.0..9.0, 0.0..max_val * max_val_margin)?;

  chart.configure_mesh().draw()?;

  // Plot actual result
  let actual_data: Vec<(f64, f64)> = x_grid.iter().zip(actual.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(actual_data, BLUE.stroke_width(actual_line_width)))?
    .label("Actual (ndarray-conv)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], BLUE));

  // Plot expected result
  let expected_data: Vec<(f64, f64)> = x_grid.iter().zip(expected.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(expected_data, RED.stroke_width(expected_line_width)))?
    .label("Expected")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], RED));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved convolution results plot to {output_path}");

  Ok(())
}

fn plot_error_analysis(
  x_grid: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
  output_dir: &str,
) -> eyre::Result<()> {
  // Plot configuration
  let plot_width = 800;
  let plot_height = 600;
  let plot_margin = 20;
  let x_label_area = 40;
  let y_label_area = 50;
  let line_width = 2;
  let legend_offset = 10;
  let max_val_margin = 1.1;
  let font_size = 24;

  let output_path = format!("{output_dir}/error_analysis.svg");
  let root = SVGBackend::new(&output_path, (plot_width, plot_height)).into_drawing_area();
  root.fill(&WHITE)?;

  // Calculate absolute error
  let absolute_errors: Vec<f64> = actual
    .iter()
    .zip(expected.iter())
    .map(|(&act, &exp)| (act - exp).abs())
    .collect();

  let max_error = absolute_errors.iter().fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption("Absolute Error: |Actual - Expected|", ("Arial", font_size))
    .margin(plot_margin)
    .x_label_area_size(x_label_area)
    .y_label_area_size(y_label_area)
    .build_cartesian_2d(-7.0..9.0, 0.0..max_error * max_val_margin)?;

  chart.configure_mesh().draw()?;

  // Plot error
  let error_data: Vec<(f64, f64)> = x_grid
    .iter()
    .zip(absolute_errors.iter())
    .map(|(&x, y)| (x, *y))
    .collect();
  chart
    .draw_series(LineSeries::new(error_data, GREEN.stroke_width(line_width)))?
    .label("Absolute Error")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], GREEN));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved error analysis plot to {output_path}");

  Ok(())
}
