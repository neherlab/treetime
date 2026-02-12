use std::cmp::min;

use eyre::Report;
use ndarray::Array1;
use num::Float;
use treetime_utils::make_error;

/// Performs piecewise linear interpolation on non-uniformly spaced data
///
/// Interpolates values from (x, y) pairs at query points using linear interpolation.
/// Produces uniformly spaced output values regardless of input spacing.
/// For query points outside the range of x, uses linear extrapolation based on the
/// nearest two points.
///
/// # Arguments
///
/// * `x` - Non-uniform x coordinates (must be sorted ascending, at least 2 points)
/// * `y` - Corresponding y values (same length as x)
/// * `n_query` - Number of query points
/// * `x_query_fn` - Closure that generates query point for given index
///
/// # Returns
///
/// Interpolated y values at query points
pub fn interp_nonuniform<T, F>(x: &Array1<T>, y: &Array1<T>, n_query: usize, x_query_fn: F) -> Result<Array1<T>, Report>
where
  T: Float,
  F: Fn(usize) -> T,
{
  if x.len() < 2 {
    return make_error!("Input arrays must have at least 2 points");
  }
  if x.len() != y.len() {
    return make_error!("x and y arrays must have the same length");
  }

  let n = x.len();
  let mut result = Array1::zeros(n_query);

  for i in 0..n_query {
    let xi = x_query_fn(i);
    if xi <= x[0] {
      let slope = (y[1] - y[0]) / (x[1] - x[0]);
      result[i] = y[0] + slope * (xi - x[0]);
    } else if xi >= x[n - 1] {
      let slope = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]);
      result[i] = y[n - 1] + slope * (xi - x[n - 1]);
    } else {
      let left = find_interval(x, xi);
      let t = (xi - x[left]) / (x[left + 1] - x[left]);
      result[i] = y[left] + t * (y[left + 1] - y[left]);
    }
  }

  Ok(result)
}

/// Finds the interval index in a sorted array containing the query value
///
/// Uses binary search to find index i such that x[i] <= xi < x[i+1].
///
/// # Arguments
///
/// * `x` - Sorted array (ascending order, at least 2 points)
/// * `xi` - Query value
///
/// # Returns
///
/// Index i where x[i] <= xi < x[i+1]. For xi >= x[n-1], returns n-2.
fn find_interval<T>(x: &Array1<T>, xi: T) -> usize
where
  T: Float,
{
  debug_assert!(x.len() >= 2, "x must have at least 2 points");
  debug_assert!(xi >= x[0] && xi <= x[x.len() - 1], "xi must be within x range");

  let n = x.len();
  let idx = x.as_slice().unwrap().partition_point(|&val| val <= xi);
  min(idx.saturating_sub(1), n - 2)
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_interp_nonuniform_exact_points() -> Result<(), Report> {
    let x = array![1.0, 2.0, 4.0, 7.0];
    let y = array![2.0, 4.0, 6.0, 8.0];
    let x_query = array![1.0, 2.0, 4.0, 7.0];

    let result = interp_nonuniform(&x, &y, x_query.len(), |i| x_query[i])?;
    assert_ulps_eq!(result, y, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_interp_nonuniform_interior() -> Result<(), Report> {
    let x = array![0.0, 1.0, 3.0];
    let y = array![0.0, 2.0, 6.0];
    let x_query = array![0.5, 2.0];

    let result = interp_nonuniform(&x, &y, x_query.len(), |i| x_query[i])?;
    let expected = array![1.0, 4.0];
    assert_ulps_eq!(result, expected, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_interp_nonuniform_left_extrapolation() -> Result<(), Report> {
    let x = array![1.0, 2.0, 3.0];
    let y = array![2.0, 4.0, 6.0];
    let x_query = array![0.0];

    let result = interp_nonuniform(&x, &y, x_query.len(), |i| x_query[i])?;
    let expected = array![0.0];
    assert_ulps_eq!(result, expected, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_interp_nonuniform_right_extrapolation() -> Result<(), Report> {
    let x = array![1.0, 2.0, 3.0];
    let y = array![2.0, 4.0, 6.0];
    let x_query = array![4.0];

    let result = interp_nonuniform(&x, &y, x_query.len(), |i| x_query[i])?;
    let expected = array![8.0];
    assert_ulps_eq!(result, expected, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_interp_nonuniform_mixed() -> Result<(), Report> {
    let x = array![1.0, 3.0, 5.0];
    let y = array![10.0, 20.0, 30.0];
    let x_query = array![0.0, 2.0, 4.0, 6.0];

    let result = interp_nonuniform(&x, &y, x_query.len(), |i| x_query[i])?;
    let expected = array![5.0, 15.0, 25.0, 35.0];
    assert_ulps_eq!(result, expected, max_ulps = 4);
    Ok(())
  }
}
