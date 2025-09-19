use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use plotters::prelude::*;
use std::f64::consts::PI;
use treetime::distribution::reference::domain_agreement_metrics::DomainAgreementMetrics;
use treetime::utils::float_fmt::float_to_digits;

fn main() -> eyre::Result<()> {
  println!("=== Gaussian Convolution Demo with ndarray-conv ===\n");

  // Configuration constants
  let output_dir = "tmp/conv_plots";

  // Parameters from Python example
  let sigma_f = 1.0_f64;
  let sigma_g = 1.0_f64;
  let mu = 1.0_f64;

  // Grid configuration
  let x_min = -7.0;
  let x_max = 9.0;
  let n_points = 321;

  // Comparison configuration
  let peak_comparison_range = 5; // ±5 points around peak

  // Create output directory
  std::fs::create_dir_all(output_dir)?;

  println!("Parameters: σ_f={sigma_f}, σ_g={sigma_g}, μ={mu}");

  // Use a narrower grid centered around the effective convolution domain
  // This allows ndarray-conv to work directly without complex interpolation
  let x_grid = Array1::linspace(x_min, x_max, n_points);
  let dx = (x_max - x_min) / (n_points - 1) as f64;

  println!("Grid: [{x_min}, {x_max}] with {n_points} points, dx={dx:.6}");

  // Define Gaussian functions on this grid
  let f_vals: Array1<f64> =
    x_grid.mapv(|x| (-0.5_f64 * (x / sigma_f).powi(2)).exp() / (sigma_f * (2.0_f64 * PI).sqrt()));

  let g_vals: Array1<f64> =
    x_grid.mapv(|x| (-0.5_f64 * ((x - mu) / sigma_g).powi(2)).exp() / (sigma_g * (2.0_f64 * PI).sqrt()));

  // Perform convolution using ndarray-conv with Full mode
  let raw_conv = f_vals.conv(&g_vals, ConvMode::Full, PaddingMode::Zeros)?;

  // Scale by grid spacing to approximate continuous convolution integral
  let conv_scaled = &raw_conv * dx;

  // The Full convolution result spans a wider domain
  // Create the extended domain for the full convolution result
  let conv_len = raw_conv.len();
  let conv_x_min = 2.0 * x_min;
  let conv_x_max = 2.0 * x_max;
  let conv_x_grid = Array1::linspace(conv_x_min, conv_x_max, conv_len);

  // Extract the portion that corresponds to our original evaluation domain
  // Find indices in conv_x_grid that correspond to our desired x_grid range
  let mut start_idx = 0;
  for (i, &conv_x) in conv_x_grid.iter().enumerate() {
    if conv_x >= x_min {
      start_idx = i;
      break;
    }
  }

  let actual = conv_scaled
    .slice(ndarray::s![start_idx..start_idx + n_points])
    .to_owned();

  // Compute expected result on the same grid
  let expected: Array1<f64> = x_grid.mapv(|x| {
    let var_sum = sigma_g.powi(2) + sigma_f.powi(2);
    (-0.5_f64 * (x - mu).powi(2) / var_sum).exp() / (var_sum * 2.0_f64 * PI).sqrt()
  });

  println!("\nndarray-conv convolution completed");
  println!("Input functions length: {}", x_grid.len());
  println!("Raw convolution length: {}", raw_conv.len());
  println!("Final result length: {}", actual.len());

  // Comprehensive domain-wide agreement metrics
  compute_domain_agreement_metrics(&x_grid, &actual, &expected)?;

  // Generate plots
  plot_input_functions(&x_grid, &f_vals, &g_vals, sigma_f, sigma_g, mu)?;
  plot_convolution_results(&x_grid, &actual, &expected)?;
  plot_error_analysis(&x_grid, &actual, &expected)?;

  // Compare values around peak
  println!("\nComparison around peak (μ = {mu}):");
  println!(
    "{:>8} {:>12} {:>12} {:>12} {:>10}",
    "x", "Actual", "Expected", "Difference", "Rel Error"
  );
  println!("{}", "-".repeat(66));

  // Find peak region on x_grid
  let peak_x = mu;
  let mut peak_idx = 0;
  for (i, &x) in x_grid.iter().enumerate() {
    if (x - peak_x).abs() < (x_grid[peak_idx] - peak_x).abs() {
      peak_idx = i;
    }
  }

  for i in (peak_idx.saturating_sub(peak_comparison_range))..(peak_idx + peak_comparison_range + 1).min(x_grid.len()) {
    let x_val = x_grid[i];
    let actual_val = actual[i];
    let expected_val = expected[i];
    let diff = (actual_val - expected_val).abs();
    let rel_err = if expected_val != 0.0 {
      diff / expected_val * 100.0
    } else {
      0.0
    };
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

  Ok(())
}

fn compute_domain_agreement_metrics(
  x_grid: &Array1<f64>,
  actual: &Array1<f64>,
  expected: &Array1<f64>,
) -> eyre::Result<()> {
  let metrics = DomainAgreementMetrics::new(x_grid, actual, expected)?;
  println!("\n{metrics}");
  Ok(())
}

fn plot_input_functions(
  x_grid: &Array1<f64>,
  f_vals: &Array1<f64>,
  g_vals: &Array1<f64>,
  sigma_f: f64,
  sigma_g: f64,
  mu: f64,
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

  let root = SVGBackend::new("tmp/conv_plots/input_functions.svg", (plot_width, plot_height)).into_drawing_area();
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
  println!("Saved input functions plot to tmp/conv_plots/input_functions.svg");

  Ok(())
}

fn plot_convolution_results(x_grid: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<()> {
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

  let root = SVGBackend::new("tmp/conv_plots/convolution_results.svg", (plot_width, plot_height)).into_drawing_area();
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
  println!("Saved convolution results plot to tmp/conv_plots/convolution_results.svg");

  Ok(())
}

fn plot_error_analysis(x_grid: &Array1<f64>, actual: &Array1<f64>, expected: &Array1<f64>) -> eyre::Result<()> {
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

  let root = SVGBackend::new("tmp/conv_plots/error_analysis.svg", (plot_width, plot_height)).into_drawing_area();
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
  println!("Saved error analysis plot to tmp/conv_plots/error_analysis.svg");

  Ok(())
}
