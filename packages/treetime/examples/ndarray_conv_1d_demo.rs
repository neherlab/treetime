use ndarray::Array1;
use ndarray_conv::{ConvExt, ConvMode, PaddingMode};
use plotters::prelude::*;
use std::f64::consts::PI;

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

  let numerical = conv_scaled
    .slice(ndarray::s![start_idx..start_idx + n_points])
    .to_owned();

  // Compute analytical result on the same grid
  let analytical: Array1<f64> = x_grid.mapv(|x| {
    let var_sum = sigma_g.powi(2) + sigma_f.powi(2);
    (-0.5_f64 * (x - mu).powi(2) / var_sum).exp() / (var_sum * 2.0_f64 * PI).sqrt()
  });

  println!("\nndarray-conv convolution completed");
  println!("Input functions length: {}", x_grid.len());
  println!("Raw convolution length: {}", raw_conv.len());
  println!("Final result length: {}", numerical.len());

  // Comprehensive domain-wide agreement metrics
  compute_domain_agreement_metrics(&x_grid, &numerical, &analytical)?;

  // Generate plots
  plot_input_functions(&x_grid, &f_vals, &g_vals, sigma_f, sigma_g, mu)?;
  plot_convolution_results(&x_grid, &numerical, &analytical)?;
  plot_error_analysis(&x_grid, &numerical, &analytical)?;

  // Compare values around peak
  println!("\nComparison around peak (μ = {mu}):");
  println!("x\t\tNumerical\tAnalytical\tDifference\tRel Error");

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
    let num_val = numerical[i];
    let ana_val = analytical[i];
    let diff = (num_val - ana_val).abs();
    let rel_err = if ana_val != 0.0 { diff / ana_val * 100.0 } else { 0.0 };
    println!("{x_val:.2}\t\t{num_val:.6}\t\t{ana_val:.6}\t\t{diff:.2e}\t\t{rel_err:.2}%");
  }

  Ok(())
}

fn compute_domain_agreement_metrics(
  x_grid: &Array1<f64>,
  numerical: &Array1<f64>,
  analytical: &Array1<f64>,
) -> eyre::Result<()> {
  // Tolerance thresholds
  let abs_tolerance_1 = 1e-6;
  let abs_tolerance_2 = 1e-9;
  let abs_tolerance_3 = 1e-12;

  let rel_tolerance_1 = 0.01; // 1%
  let rel_tolerance_2 = 0.001; // 0.1%
  let rel_tolerance_3 = 0.0001; // 0.01%

  // R² quality thresholds
  let r2_excellent = 0.999999;
  let r2_very_good = 0.9999;
  let r2_good = 0.99;

  println!("\n=== Domain-Wide Agreement Metrics ===");

  // Basic error statistics
  let abs_errors: Vec<f64> = numerical
    .iter()
    .zip(analytical.iter())
    .map(|(&num, &ana)| (num - ana).abs())
    .collect();

  let rel_errors: Vec<f64> = numerical
    .iter()
    .zip(analytical.iter())
    .map(|(&num, &ana)| {
      if ana.abs() > 1e-15 {
        (num - ana).abs() / ana.abs()
      } else {
        0.0
      }
    })
    .collect();

  let mean_abs_error = abs_errors.iter().sum::<f64>() / abs_errors.len() as f64;
  let max_abs_error = abs_errors.iter().fold(0.0_f64, |a, &b| a.max(b));
  let std_abs_error = {
    let variance = abs_errors.iter().map(|&x| (x - mean_abs_error).powi(2)).sum::<f64>() / abs_errors.len() as f64;
    variance.sqrt()
  };

  let mean_rel_error = rel_errors.iter().sum::<f64>() / rel_errors.len() as f64;
  let max_rel_error = rel_errors.iter().fold(0.0_f64, |a, &b| a.max(b));

  // R-squared (coefficient of determination)
  let mean_analytical = analytical.iter().sum::<f64>() / analytical.len() as f64;
  let ss_tot: f64 = analytical.iter().map(|&ana| (ana - mean_analytical).powi(2)).sum();
  let ss_res: f64 = numerical
    .iter()
    .zip(analytical.iter())
    .map(|(&num, &ana)| (ana - num).powi(2))
    .sum();
  let r_squared = 1.0 - (ss_res / ss_tot);

  // Correlation coefficient
  let mean_numerical = numerical.iter().sum::<f64>() / numerical.len() as f64;

  let numerator: f64 = numerical
    .iter()
    .zip(analytical.iter())
    .map(|(&num, &ana)| (num - mean_numerical) * (ana - mean_analytical))
    .sum();

  let denom_num = numerical
    .iter()
    .map(|&num| (num - mean_numerical).powi(2))
    .sum::<f64>()
    .sqrt();

  let denom_ana = analytical
    .iter()
    .map(|&ana| (ana - mean_analytical).powi(2))
    .sum::<f64>()
    .sqrt();

  let correlation = numerator / (denom_num * denom_ana);

  // Root Mean Square Error (RMSE)
  let rmse = (abs_errors.iter().map(|&x| x.powi(2)).sum::<f64>() / abs_errors.len() as f64).sqrt();

  // Mean Absolute Percentage Error (MAPE)
  let mape = rel_errors.iter().sum::<f64>() / rel_errors.len() as f64 * 100.0;

  // Points within tolerance thresholds
  let within_1e_6 = abs_errors.iter().filter(|&&x| x < abs_tolerance_1).count();
  let within_1e_9 = abs_errors.iter().filter(|&&x| x < abs_tolerance_2).count();
  let within_1e_12 = abs_errors.iter().filter(|&&x| x < abs_tolerance_3).count();

  let within_1_percent = rel_errors.iter().filter(|&&x| x < rel_tolerance_1).count();
  let within_0_1_percent = rel_errors.iter().filter(|&&x| x < rel_tolerance_2).count();
  let within_0_01_percent = rel_errors.iter().filter(|&&x| x < rel_tolerance_3).count();

  let total_points = abs_errors.len();

  // Peak error location
  let max_error_idx = abs_errors
    .iter()
    .enumerate()
    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
    .map_or(0, |(i, _)| i);

  // Print comprehensive metrics
  println!("Total evaluation points: {total_points}");
  println!();

  println!("Absolute Error Statistics:");
  println!("  Mean absolute error:    {mean_abs_error:.6e}");
  println!("  Max absolute error:     {max_abs_error:.6e}");
  println!("  Std absolute error:     {std_abs_error:.6e}");
  println!("  RMSE:                   {rmse:.6e}");
  println!("  Max error at x = {:.3}", x_grid[max_error_idx]);
  println!();

  println!("Relative Error Statistics:");
  println!(
    "  Mean relative error:    {:.6e} ({:.4}%)",
    mean_rel_error,
    mean_rel_error * 100.0
  );
  println!(
    "  Max relative error:     {:.6e} ({:.4}%)",
    max_rel_error,
    max_rel_error * 100.0
  );
  println!("  MAPE:                   {mape:.6}%");
  println!();

  println!("Agreement Quality:");
  println!("  R² (coefficient of determination): {r_squared:.12}");
  println!("  Correlation coefficient:            {correlation:.12}");
  println!();

  println!("Points within absolute tolerance:");
  println!(
    "  < {:.0e}:  {:3}/{} ({:5.1}%)",
    abs_tolerance_1,
    within_1e_6,
    total_points,
    within_1e_6 as f64 / total_points as f64 * 100.0
  );
  println!(
    "  < {:.0e}:  {:3}/{} ({:5.1}%)",
    abs_tolerance_2,
    within_1e_9,
    total_points,
    within_1e_9 as f64 / total_points as f64 * 100.0
  );
  println!(
    "  < {:.0e}: {:3}/{} ({:5.1}%)",
    abs_tolerance_3,
    within_1e_12,
    total_points,
    within_1e_12 as f64 / total_points as f64 * 100.0
  );
  println!();

  println!("Points within relative tolerance:");
  println!(
    "  < {:.0}%:    {:3}/{} ({:5.1}%)",
    rel_tolerance_1 * 100.0,
    within_1_percent,
    total_points,
    within_1_percent as f64 / total_points as f64 * 100.0
  );
  println!(
    "  < {:.1}%:  {:3}/{} ({:5.1}%)",
    rel_tolerance_2 * 100.0,
    within_0_1_percent,
    total_points,
    within_0_1_percent as f64 / total_points as f64 * 100.0
  );
  println!(
    "  < {:.2}%: {:3}/{} ({:5.1}%)",
    rel_tolerance_3 * 100.0,
    within_0_01_percent,
    total_points,
    within_0_01_percent as f64 / total_points as f64 * 100.0
  );

  // Overall assessment
  println!();
  if r_squared > r2_excellent {
    println!("🟢 EXCELLENT: Near-perfect agreement (R² > {r2_excellent:.6})");
  } else if r_squared > r2_very_good {
    println!("🟡 VERY GOOD: High agreement (R² > {r2_very_good:.4})");
  } else if r_squared > r2_good {
    println!("🟠 GOOD: Reasonable agreement (R² > {r2_good:.2})");
  } else {
    println!("🔴 POOR: Low agreement (R² < {r2_good:.2})");
  }

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

fn plot_convolution_results(
  x_grid: &Array1<f64>,
  numerical: &Array1<f64>,
  analytical: &Array1<f64>,
) -> eyre::Result<()> {
  // Plot configuration
  let plot_width = 800;
  let plot_height = 600;
  let plot_margin = 20;
  let x_label_area = 40;
  let y_label_area = 50;
  let numerical_line_width = 2;
  let analytical_line_width = 3;
  let legend_offset = 10;
  let max_val_margin = 1.1;
  let font_size = 24;

  let root = SVGBackend::new("tmp/conv_plots/convolution_results.svg", (plot_width, plot_height)).into_drawing_area();
  root.fill(&WHITE)?;

  let max_val = numerical
    .iter()
    .chain(analytical.iter())
    .fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption("Convolution Results: (f * g)(x)", ("Arial", font_size))
    .margin(plot_margin)
    .x_label_area_size(x_label_area)
    .y_label_area_size(y_label_area)
    .build_cartesian_2d(-7.0..9.0, 0.0..max_val * max_val_margin)?;

  chart.configure_mesh().draw()?;

  // Plot numerical result
  let num_data: Vec<(f64, f64)> = x_grid.iter().zip(numerical.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(num_data, BLUE.stroke_width(numerical_line_width)))?
    .label("Numerical (ndarray-conv)")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], BLUE));

  // Plot analytical result
  let ana_data: Vec<(f64, f64)> = x_grid.iter().zip(analytical.iter()).map(|(&x, &y)| (x, y)).collect();
  chart
    .draw_series(LineSeries::new(ana_data, RED.stroke_width(analytical_line_width)))?
    .label("Analytical Solution")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], RED));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved convolution results plot to tmp/conv_plots/convolution_results.svg");

  Ok(())
}

fn plot_error_analysis(x_grid: &Array1<f64>, numerical: &Array1<f64>, analytical: &Array1<f64>) -> eyre::Result<()> {
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
  let abs_error: Vec<f64> = numerical
    .iter()
    .zip(analytical.iter())
    .map(|(&num, &ana)| (num - ana).abs())
    .collect();

  let max_error = abs_error.iter().fold(0.0_f64, |a, &b| a.max(b));

  let mut chart = ChartBuilder::on(&root)
    .caption("Absolute Error: |Numerical - Analytical|", ("Arial", font_size))
    .margin(plot_margin)
    .x_label_area_size(x_label_area)
    .y_label_area_size(y_label_area)
    .build_cartesian_2d(-7.0..9.0, 0.0..max_error * max_val_margin)?;

  chart.configure_mesh().draw()?;

  // Plot error
  let error_data: Vec<(f64, f64)> = x_grid.iter().zip(abs_error.iter()).map(|(&x, y)| (x, *y)).collect();
  chart
    .draw_series(LineSeries::new(error_data, GREEN.stroke_width(line_width)))?
    .label("Absolute Error")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + legend_offset, y)], GREEN));

  chart.configure_series_labels().draw()?;
  root.present()?;
  println!("Saved error analysis plot to tmp/conv_plots/error_analysis.svg");

  Ok(())
}
