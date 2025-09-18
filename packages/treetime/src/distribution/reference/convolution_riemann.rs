use crate::distribution::reference::grid_fn::GridFn;
use eyre::Report;
use ndarray::Array1;

/// Compute convolution using Riemann sum approach (exactly like Python version)
///
/// Computes (f * g)(x) = ∫ f(s) * g(x-s) ds using discrete Riemann sum:
/// ds * Σ f(s_i) * g(x - s_i)
///
/// # Arguments
/// * `f` - First function
/// * `g` - Second function
/// * `x` - Evaluation point
/// * `s_grid` - Integration grid (must be uniformly spaced)
///
/// # Returns
/// Convolution value at point x
pub fn convolve_riemann(f: &GridFn, g: &GridFn, x: f64, s_grid: &Array1<f64>) -> Result<f64, Report> {
  assert!(s_grid.len() >= 2, "Integration grid must have at least 2 points");

  // Calculate grid spacing (assume uniform grid)
  let ds = s_grid[1] - s_grid[0];

  // Evaluate f(s_grid) and g(x - s_grid)
  let mut f_vals = Array1::zeros(s_grid.len());
  let mut g_vals = Array1::zeros(s_grid.len());

  for (i, &s) in s_grid.iter().enumerate() {
    f_vals[i] = f.interp(s).unwrap_or(0.0); // Use 0 for out-of-domain
    g_vals[i] = g.interp(x - s).unwrap_or(0.0); // Use 0 for out-of-domain
  }

  // Compute Riemann sum: ds * sum(f(s_grid) * g(x - s_grid))
  Ok(ds * (&f_vals * &g_vals).sum())
}

/// Compute convolution using shared grids (no interpolation needed)
/// This is more efficient when f, g, and result all use the same grid
pub fn convolve_riemann_shared_grid(
  f: &GridFn,
  g: &GridFn,
  x_grid: &Array1<f64>,
  s_grid: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  assert!(s_grid.len() >= 2, "Integration grid must have at least 2 points");

  // Calculate grid spacing (assume uniform grid)
  let ds = s_grid[1] - s_grid[0];

  // Pre-evaluate f on s_grid (avoid interpolation if possible)
  let f_on_s_grid = f.eval_array(s_grid).unwrap_or_else(|_| Array1::zeros(s_grid.len()));

  let mut result = Array1::zeros(x_grid.len());

  for (i, &x) in x_grid.iter().enumerate() {
    // Create x - s_grid array
    let x_minus_s: Array1<f64> = s_grid.mapv(|s| x - s);

    // Evaluate g at x - s points
    let g_vals = g.eval_array(&x_minus_s).unwrap_or_else(|_| Array1::zeros(s_grid.len()));

    // Compute Riemann sum: ds * sum(f(s_grid) * g(x - s_grid))
    result[i] = ds * (&f_on_s_grid * &g_vals).sum();
  }

  Ok(result)
}

/// Compute convolution over multiple x values
pub fn convolve_riemann_grid(
  f: &GridFn,
  g: &GridFn,
  x_grid: &Array1<f64>,
  s_grid: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  let mut result = Array1::zeros(x_grid.len());
  for (i, &x) in x_grid.iter().enumerate() {
    result[i] = convolve_riemann(f, g, x, s_grid)?;
  }
  Ok(result)
}

/// Create GridFn from convolution result
pub fn convolve_riemann_as_gridfn(
  f: &GridFn,
  g: &GridFn,
  x_grid: Array1<f64>,
  s_grid: &Array1<f64>,
) -> Result<GridFn, Report> {
  let y_values = convolve_riemann_grid(f, g, &x_grid, s_grid)?;
  GridFn::new(x_grid, y_values)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::distribution::reference::exponential::{
    exponential_convolution_on_grid, exponential_f_on_grid, exponential_g_on_grid,
  };
  use crate::distribution::reference::gaussian::{gaussian_convolution, gaussian_f, gaussian_g};
  use crate::pretty_assert_ulps_eq;
  use approx::assert_ulps_eq;
  use ndarray::Array1;

  #[test]
  fn test_gaussian_convolution() -> Result<(), Report> {
    // Test convolution of two Gaussians
    let sigma_f = 1.0;
    let sigma_g = 1.0;
    let mu = 1.0;

    // Use wider domains to avoid interpolation issues during convolution
    let f = gaussian_f(sigma_f, (-10.0, 10.0), 0.1)?;
    let g = gaussian_g(sigma_g, mu, (-10.0, 10.0), 0.1)?;

    // Integration grid covering full domain
    let s_grid = Array1::from_iter((0..201).map(|i| -10.0 + i as f64 * 0.1));

    // Test single point convolution
    let x = 0.0;
    let result = convolve_riemann(&f, &g, x, &s_grid)?;

    // Use existing analytical convolution function
    let analytical_fn = gaussian_convolution(sigma_f, sigma_g, mu, (-1.0, 1.0), 0.1)?;
    let expected = analytical_fn.interp(x)?;

    assert_ulps_eq!(result, expected, max_ulps = 1000);
    Ok(())
  }

  #[test]
  fn test_exponential_convolution() -> Result<(), Report> {
    // Test convolution of two exponentials
    let a_f = 1.0;
    let a_g = 2.0;

    // Create shared grids to avoid interpolation
    let s_grid = Array1::from_iter((0..2001).map(|i| -5.0 + i as f64 * 0.01));
    let x_grid = Array1::from_iter((0..51).map(|i| -1.0 + i as f64 * 0.1));

    // Create functions on shared grids
    let f = exponential_f_on_grid(a_f, &s_grid)?;
    let g = exponential_g_on_grid(a_g, &s_grid)?;

    // Compute convolution using shared grid approach
    let result = convolve_riemann_shared_grid(&f, &g, &x_grid, &s_grid)?;

    // Get analytical solution on the same x_grid
    let analytical_fn = exponential_convolution_on_grid(a_f, a_g, &x_grid)?;
    let expected = analytical_fn.y();

    // Compare all points using pretty_assert_ulps_eq! (no interpolation needed)
    pretty_assert_ulps_eq!(&result, expected, epsilon = 1e-1);
    Ok(())
  }

  #[test]
  fn test_convolution_grid() -> Result<(), Report> {
    // Test grid convolution functionality with wider domains
    let f = gaussian_f(1.0, (-6.0, 6.0), 0.1)?;
    let g = gaussian_f(1.0, (-6.0, 6.0), 0.1)?;

    let x_grid = Array1::from_iter((0..11).map(|i| -1.0 + i as f64 * 0.2));
    let s_grid = Array1::from_iter((0..121).map(|i| -6.0 + i as f64 * 0.1));

    let result = convolve_riemann_grid(&f, &g, &x_grid, &s_grid)?;

    assert_eq!(result.len(), x_grid.len());

    // Result should be symmetric around 0 for identical centered Gaussians
    let mid_idx = result.len() / 2;
    let max_val = result[mid_idx];
    assert!(max_val > 0.0);

    Ok(())
  }

  #[test]
  fn test_convolution_as_gridfn() -> Result<(), Report> {
    // Test creating GridFn from convolution with wider domains
    let f = gaussian_f(1.0, (-5.0, 5.0), 0.1)?;
    let g = gaussian_f(1.0, (-5.0, 5.0), 0.1)?;

    let x_grid = Array1::from_iter((0..21).map(|i| -2.0 + i as f64 * 0.2));
    let s_grid = Array1::from_iter((0..101).map(|i| -5.0 + i as f64 * 0.1));

    let conv_fn = convolve_riemann_as_gridfn(&f, &g, x_grid.clone(), &s_grid)?;

    pretty_assert_ulps_eq!(conv_fn.x(), &x_grid, max_ulps = 1);
    assert_eq!(conv_fn.y().len(), x_grid.len());

    // Test interpolation works
    let interp_val = conv_fn.interp(0.0)?;
    assert!(interp_val > 0.0);

    Ok(())
  }
}
