use eyre::Report;
use ndarray::Array1;

use treetime::distribution::reference::convolution_riemann::{
  convolve_riemann, convolve_riemann_as_gridfn, convolve_riemann_grid,
};
use treetime::distribution::reference::exponential::{exponential_convolution, exponential_f, exponential_g};
use treetime::distribution::reference::gaussian::{gaussian_convolution, gaussian_f, gaussian_g};

fn main() -> Result<(), Report> {
  println!("🧮 TreeTime Convolution Demo using Riemann Sum Approach");
  println!("{}", "=".repeat(60));

  // Example 1: Gaussian Convolution
  println!("\n📊 Example 1: Convolution of Two Gaussians");
  println!("{}", "-".repeat(40));

  let sigma_f = 1.0;
  let sigma_g = 1.5;
  let mu = 2.0;

  println!("Parameters:");
  println!("  f(x) = Gaussian(σ={sigma_f}, μ=0)");
  println!("  g(x) = Gaussian(σ={sigma_g}, μ={mu})");

  // Create Gaussian functions with wide domains
  let f = gaussian_f(sigma_f, (-8.0, 8.0), 0.05)?;
  let g = gaussian_g(sigma_g, mu, (-8.0, 8.0), 0.05)?;

  // Integration grid
  let s_grid = Array1::from_iter((0..321).map(|i| -8.0 + i as f64 * 0.05));

  // Test single point convolution
  let test_points = vec![0.0, 1.0, 2.0, 3.0];
  // Create analytical convolution function for comparison
  let analytical_fn = gaussian_convolution(sigma_f, sigma_g, mu, (-10.0, 10.0), 0.05)?;

  println!("\nSingle point convolutions:");
  for x in test_points {
    let numerical = convolve_riemann(&f, &g, x, &s_grid)?;
    let analytical = analytical_fn.interp(x)?;
    let error = ((numerical - analytical) / analytical * 100.0).abs();

    println!("  x={x:.1}: numerical={numerical:.6}, analytical={analytical:.6}, error={error:.3}%");
  }

  // Example 2: Exponential Convolution
  println!("\n📊 Example 2: Convolution of Two Exponentials");
  println!("{}", "-".repeat(40));

  let a_f = 1.0;
  let a_g = 2.0;

  println!("Parameters:");
  println!("  f(x) = {a_f} * exp(-{a_f}x) for x≥0");
  println!("  g(x) = {a_g} * exp(-{a_g}x) for x≥0");

  // Create exponential functions
  let f_exp = exponential_f(a_f, (-2.0, 10.0), 0.01)?;
  let g_exp = exponential_g(a_g, (-2.0, 10.0), 0.01)?;

  // Integration grid for exponentials
  let s_grid_exp = Array1::from_iter((0..1201).map(|i| -2.0 + i as f64 * 0.01));

  // Create analytical convolution function for comparison
  let analytical_exp_fn = exponential_convolution(a_f, a_g, (-2.0, 10.0), 0.01)?;

  println!("\nSingle point convolutions:");
  let test_points_exp = vec![0.5, 1.0, 2.0, 3.0];
  for x in test_points_exp {
    let numerical = convolve_riemann(&f_exp, &g_exp, x, &s_grid_exp)?;
    let analytical = analytical_exp_fn.interp(x)?;
    let error = ((numerical - analytical) / analytical * 100.0).abs();

    println!("  x={x:.1}: numerical={numerical:.6}, analytical={analytical:.6}, error={error:.3}%");
  }

  // Example 3: Grid Convolution
  println!("\n📊 Example 3: Grid-based Convolution");
  println!("{}", "-".repeat(40));

  // Evaluate convolution over a grid of points
  let x_grid = Array1::from_iter((0..21).map(|i| -2.0 + i as f64 * 0.2));
  let result_grid = convolve_riemann_grid(&f, &g, &x_grid, &s_grid)?;

  println!("Convolution evaluated on {} grid points:", x_grid.len());
  for (i, (&x, &y)) in x_grid.iter().zip(result_grid.iter()).enumerate() {
    if i % 5 == 0 {
      // Print every 5th point to avoid clutter
      println!("  x={x:.1}: y={y:.6}");
    }
  }

  // Example 4: Creating GridFn from Convolution
  println!("\n📊 Example 4: GridFn Convolution Result");
  println!("{}", "-".repeat(40));

  let conv_fn = convolve_riemann_as_gridfn(&f, &g, x_grid, &s_grid)?;

  println!("Created GridFn with {} points", conv_fn.x().len());
  println!(
    "Maximum value: {:.6} at x={:.2}",
    conv_fn.max_value(),
    conv_fn.max_position()
  );

  // Test interpolation
  let interp_test_x = 1.5;
  let interp_val = conv_fn.interp(interp_test_x)?;
  println!("Interpolated value at x={interp_test_x}: {interp_val:.6}");

  // Performance comparison
  println!("\n⚡ Performance Note:");
  println!("The Riemann sum approach provides:");
  println!("  ✓ Simple implementation");
  println!("  ✓ Easy to understand and debug");
  println!("  ✓ Works with any function pair");
  println!("  ⚠ Lower accuracy than analytical solutions");
  println!("  ⚠ Computational cost scales with grid resolution");

  println!("\n✅ Demo completed successfully!");

  Ok(())
}
