use std::cmp::min;

use crate::InterpElem;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_utils::serde::{array1_as_vec, array1_from_vec};

/// Function represented on a regular grid for piecewise linear interpolation
///
/// Represents a function as a set of (x, y) points on a grid, providing linear interpolation
/// between points and linear extrapolation beyond the grid boundaries.
///
/// # Invariants
///
/// - `x` array must be strictly monotonically increasing: `x[i] < x[i+1]` for all `i`
/// - `x` and `y` arrays must have equal length
/// - Both arrays must contain at least 2 points
/// - These invariants are checked via `debug_assert!` in the constructor
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound(serialize = "T: Serialize", deserialize = "T: Deserialize<'de>"))]
pub struct GridFn<T: InterpElem> {
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  x: Array1<T>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  y: Array1<T>,
}

impl<T: InterpElem> GridFn<T> {
  pub fn new(x: Array1<T>, y: Array1<T>) -> Result<Self, Report> {
    debug_assert_eq!(x.len(), y.len());
    debug_assert!(x.len() >= 2);
    debug_assert!(
      x.windows(2).into_iter().all(|w| w[0] < w[1]),
      "x array must be strictly monotonic increasing"
    );

    Ok(Self { x, y })
  }

  pub fn from_grid<F>((x_min, x_max): (T, T), dx: T, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    debug_assert!(dx > T::zero());
    debug_assert!(x_max > x_min);

    let n_points = ((x_max - x_min) / dx + T::one()).round().to_usize().unwrap();
    Self::from_n_points((x_min, x_max), n_points, y_fn)
  }

  pub fn from_n_points<F>((x_min, x_max): (T, T), n_points: usize, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    debug_assert!(n_points >= 2);
    debug_assert!(x_max > x_min);

    let x = Array1::linspace(x_min, x_max, n_points);
    let y = x.mapv(y_fn);
    Self::new(x, y)
  }

  pub fn from_fn_samples<F>(x_samples: Array1<T>, y_fn: F) -> Result<Self, Report>
  where
    F: Fn(T) -> T,
  {
    let y = x_samples.mapv(y_fn);
    Self::new(x_samples, y)
  }

  pub fn constant((x_min, x_max): (T, T), n_points: usize, value: T) -> Result<Self, Report>
  where
    T: Float,
  {
    debug_assert!(n_points >= 2);
    debug_assert!(x_max > x_min);

    let x = Array1::linspace(x_min, x_max, n_points);
    let y = Array1::from_elem(n_points, value);
    Self::new(x, y)
  }

  pub fn zeros((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::constant((x_min, x_max), n_points, T::zero())
  }

  pub fn ones((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::constant((x_min, x_max), n_points, T::one())
  }

  pub fn from_pairs(pairs: Vec<(T, T)>) -> Result<Self, Report> {
    debug_assert!(pairs.len() >= 2);

    let (x, y): (Vec<T>, Vec<T>) = pairs.into_iter().unzip();
    let x = Array1::from_vec(x);
    let y = Array1::from_vec(y);
    Self::new(x, y)
  }

  pub fn x(&self) -> &Array1<T> {
    &self.x
  }

  pub fn x_mut(&mut self) -> &mut Array1<T> {
    &mut self.x
  }

  pub fn y(&self) -> &Array1<T> {
    &self.y
  }

  pub fn y_mut(&mut self) -> &mut Array1<T> {
    &mut self.y
  }

  pub fn x_min(&self) -> T {
    self.x[0]
  }

  pub fn x_max(&self) -> T {
    self.x[self.x.len() - 1]
  }

  pub fn x_range(&self) -> (T, T) {
    (self.x_min(), self.x_max())
  }

  pub fn n_points(&self) -> usize {
    self.x.len()
  }

  pub fn dx(&self) -> T
  where
    T: Float,
  {
    (self.x_max() - self.x_min()) / T::from(self.n_points() - 1).unwrap()
  }

  pub fn y_min(&self) -> T {
    *self.y.min().unwrap()
  }

  pub fn y_max(&self) -> T {
    *self.y.max().unwrap()
  }

  pub fn y_range(&self) -> (T, T) {
    (self.y_min(), self.y_max())
  }

  pub fn to_pairs(&self) -> Vec<(T, T)> {
    self.x.iter().zip(self.y.iter()).map(|(&x, &y)| (x, y)).collect()
  }

  /// Interpolate function value at a single point
  ///
  /// Uses piecewise linear interpolation within grid bounds and linear extrapolation
  /// outside bounds using slopes from the nearest grid intervals.
  ///
  /// # Arguments
  ///
  /// * `xi` - Point at which to evaluate the function
  ///
  /// # Returns
  ///
  /// Interpolated or extrapolated function value at `xi`
  pub fn interp(&self, xi: T) -> Result<T, Report> {
    let n = self.x.len();
    let x_min = self.x[0];
    let x_max = self.x[n - 1];

    if xi < x_min {
      return Ok(self.extrapolate_left(xi));
    }

    if xi > x_max {
      return Ok(self.extrapolate_right(xi));
    }

    let idx = binary_search_interval(&self.x, xi);
    Ok(self.interpolate_at(xi, idx))
  }

  /// Interpolate function values at multiple points using optimized two-pointer algorithm
  ///
  /// Uses piecewise linear interpolation within grid bounds and linear extrapolation
  /// outside bounds. Optimized for sorted query points using a single linear scan
  /// through the grid (O(n + m) complexity where n = query points, m = grid points).
  ///
  /// # Arguments
  ///
  /// * `queries` - Query points at which to evaluate the function
  ///
  /// # Returns
  ///
  /// Array of interpolated or extrapolated function values at each query point
  ///
  /// # Preconditions
  ///
  /// * `queries` must be sorted in non-decreasing order (checked via `debug_assert!`)
  pub fn interp_many(&self, queries: &Array1<T>) -> Result<Array1<T>, Report> {
    debug_assert!(
      queries.windows(2).into_iter().all(|w| w[0] <= w[1]),
      "queries must be sorted in non-decreasing order"
    );

    let n = self.x.len();
    let x_min = self.x[0];
    let x_max = self.x[n - 1];
    let mut result = Array1::zeros(queries.len());
    let mut grid_idx = 0_usize;

    for (i, &query) in queries.iter().enumerate() {
      if query < x_min {
        result[i] = self.extrapolate_left(query);
      } else if query > x_max {
        result[i] = self.extrapolate_right(query);
      } else {
        while grid_idx < n - 1 && self.x[grid_idx + 1] < query {
          grid_idx += 1;
        }
        let idx = min(grid_idx, n - 2);
        result[i] = self.interpolate_at(query, idx);
      }
    }
    Ok(result)
  }

  fn extrapolate_left(&self, q: T) -> T {
    let slope = (self.y[1] - self.y[0]) / (self.x[1] - self.x[0]);
    self.y[0] + slope * (q - self.x[0])
  }

  fn extrapolate_right(&self, q: T) -> T {
    let n = self.x.len();
    let slope = (self.y[n - 1] - self.y[n - 2]) / (self.x[n - 1] - self.x[n - 2]);
    self.y[n - 1] + slope * (q - self.x[n - 1])
  }

  fn interpolate_at(&self, q: T, idx: usize) -> T {
    let n = self.x.len();
    if idx >= n - 1 {
      return self.y[n - 1];
    }
    let x0 = self.x[idx];
    let x1 = self.x[idx + 1];
    let y0 = self.y[idx];
    let y1 = self.y[idx + 1];
    let t = (q - x0) / (x1 - x0);
    y0 + t * (y1 - y0)
  }
}

fn binary_search_interval<T: InterpElem>(x: &Array1<T>, xi: T) -> usize {
  let n = x.len();
  if n == 0 {
    return 0;
  }
  if xi <= x[0] {
    return 0;
  }
  if xi >= x[n - 1] {
    return n - 1;
  }

  let mut left = 0_usize;
  let mut right = n - 1;

  while right - left > 1 {
    let mid = usize::midpoint(left, right);
    if xi < x[mid] {
      right = mid;
    } else {
      left = mid;
    }
  }

  left
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::{assert_abs_diff_eq, assert_ulps_eq};
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 0.0, 0.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 1.0, 10.0)]
  #[case(array![0.0, 1.0, 2.0], array![0.0, 5.0, 10.0], 1.0, 5.0)]
  #[case(array![-5.0, 5.0], array![100.0, 200.0], -5.0, 100.0)]
  #[case(array![-5.0, 5.0], array![100.0, 200.0], 5.0, 200.0)]
  #[trace]
  fn test_gridfn_interp_exact_grid_points(
    #[case] x: Array1<f64>,
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::new(x, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 0.5, 5.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 0.3, 3.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 0.7, 7.0)]
  #[case(array![0.0, 2.0], array![10.0, 20.0], 1.0, 15.0)]
  #[case(array![-1.0, 1.0], array![0.0, 100.0], 0.0, 50.0)]
  #[case(array![0.0, 1.0, 2.0], array![0.0, 10.0, 20.0], 0.5, 5.0)]
  #[case(array![0.0, 1.0, 2.0], array![0.0, 10.0, 20.0], 1.5, 15.0)]
  #[trace]
  fn test_gridfn_interp_interior_points(
    #[case] x: Array1<f64>,
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::new(x, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case(array![0.0, 1.0], array![0.0, 10.0], -1.0, -10.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], -0.5, -5.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], -2.0, -20.0)]
  #[case(array![1.0, 2.0], array![5.0, 15.0], 0.0, -5.0)]
  #[case(array![1.0, 2.0, 3.0], array![10.0, 20.0, 30.0], 0.0, 0.0)]
  #[trace]
  fn test_gridfn_interp_left_extrapolation(
    #[case] x: Array1<f64>,
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::new(x, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 2.0, 20.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 1.5, 15.0)]
  #[case(array![0.0, 1.0], array![0.0, 10.0], 3.0, 30.0)]
  #[case(array![1.0, 2.0], array![5.0, 15.0], 3.0, 25.0)]
  #[case(array![0.0, 1.0, 2.0], array![10.0, 20.0, 30.0], 3.0, 40.0)]
  #[trace]
  fn test_gridfn_interp_right_extrapolation(
    #[case] x: Array1<f64>,
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::new(x, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case(array![0.0, 1.0], array![5.0, 5.0], 0.5, 5.0)]
  #[case(array![0.0, 1.0], array![5.0, 5.0], -1.0, 5.0)]
  #[case(array![0.0, 1.0], array![5.0, 5.0], 2.0, 5.0)]
  #[case(array![0.0, 1.0, 2.0], array![42.0, 42.0, 42.0], 1.5, 42.0)]
  #[trace]
  fn test_gridfn_interp_constant_function(
    #[case] x: Array1<f64>,
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::new(x, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case(array![0.0, 1.0], array![0.0, 1e-15], 0.5, 5e-16)]
  #[case(array![0.0, 1e-10], array![0.0, 1.0], 5e-11, 0.5)]
  #[case(array![1e10, 1e10 + 1.0], array![0.0, 1.0], 1e10 + 0.5, 0.5)]
  #[trace]
  fn test_gridfn_interp_numeric_precision(
    #[case] x: Array1<f64>,
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::new(x, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0], array![0.0, 10.0, 20.0])?;
    let queries = array![0.5, 1.0, 1.5];
    let expected = array![5.0, 10.0, 15.0];
    let actual = grid_fn.interp_many(&queries)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_grid() -> Result<(), Report> {
    let grid_fn = GridFn::from_grid((0.0, 1.0), 0.25, |x| x * x)?;
    assert_eq!(grid_fn.x().len(), 5);
    assert_ulps_eq!(grid_fn.x()[0], 0.0, epsilon = 1e-14);
    assert_ulps_eq!(grid_fn.x()[4], 1.0, epsilon = 1e-14);
    assert_ulps_eq!(grid_fn.y()[0], 0.0, epsilon = 1e-14);
    assert_ulps_eq!(grid_fn.y()[2], 0.25, epsilon = 1e-14);
    assert_ulps_eq!(grid_fn.y()[4], 1.0, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_accessors() -> Result<(), Report> {
    let x = array![0.0, 1.0, 2.0];
    let y = array![10.0, 20.0, 30.0];
    let grid_fn = GridFn::new(x.clone(), y.clone())?;
    assert_eq!(&x, grid_fn.x());
    assert_eq!(&y, grid_fn.y());
    Ok(())
  }

  #[test]
  fn test_gridfn_x_min_x_max() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0, 3.0], array![10.0, 20.0, 30.0, 40.0])?;
    assert_ulps_eq!(grid_fn.x_min(), 0.0);
    assert_ulps_eq!(grid_fn.x_max(), 3.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_x_range() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![-5.0, 0.0, 5.0], array![100.0, 50.0, 0.0])?;
    let expected = (-5.0, 5.0);
    let actual = grid_fn.x_range();
    assert_ulps_eq!(expected.0, actual.0);
    assert_ulps_eq!(expected.1, actual.1);
    Ok(())
  }

  #[test]
  fn test_gridfn_n_points() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 10.0), 25, |x| x)?;
    assert_eq!(grid_fn.n_points(), 25);
    Ok(())
  }

  #[test]
  fn test_gridfn_dx_uniform() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 10.0), 11, |x| x)?;
    let expected = 1.0;
    let actual = grid_fn.dx();
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_dx_from_grid() -> Result<(), Report> {
    let dx = 0.25;
    let grid_fn = GridFn::from_grid((0.0, 1.0), dx, |x| x * x)?;
    assert_ulps_eq!(grid_fn.dx(), dx, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_parameters_roundtrip() -> Result<(), Report> {
    let x_range = (5.0, 15.0);
    let n_points = 50;
    let grid_fn = GridFn::from_n_points(x_range, n_points, |x| x.sin())?;
    assert_ulps_eq!(grid_fn.x_min(), x_range.0);
    assert_ulps_eq!(grid_fn.x_max(), x_range.1);
    assert_eq!(grid_fn.n_points(), n_points);
    assert_ulps_eq!(grid_fn.x_range().0, x_range.0);
    assert_ulps_eq!(grid_fn.x_range().1, x_range.1);
    Ok(())
  }

  #[test]
  fn test_gridfn_dx_negative_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((-10.0, -5.0), 6, |x| x)?;
    let expected = 1.0;
    let actual = grid_fn.dx();
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[test]
  fn test_gridfn_y_min_y_max() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0, 3.0], array![10.0, 5.0, 30.0, 15.0])?;
    assert_ulps_eq!(grid_fn.y_min(), 5.0);
    assert_ulps_eq!(grid_fn.y_max(), 30.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_y_range() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0, 3.0], array![-10.0, 20.0, 5.0, -5.0])?;
    let expected = (-10.0, 20.0);
    let actual = grid_fn.y_range();
    assert_ulps_eq!(expected.0, actual.0);
    assert_ulps_eq!(expected.1, actual.1);
    Ok(())
  }

  #[test]
  fn test_gridfn_y_range_constant() -> Result<(), Report> {
    let grid_fn = GridFn::constant((0.0, 10.0), 20, 42.0)?;
    assert_ulps_eq!(grid_fn.y_min(), 42.0);
    assert_ulps_eq!(grid_fn.y_max(), 42.0);
    let y_range = grid_fn.y_range();
    assert_ulps_eq!(y_range.0, 42.0);
    assert_ulps_eq!(y_range.1, 42.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_y_range_monotonic() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 10.0), 11, |x| x)?;
    assert_ulps_eq!(grid_fn.y_min(), 0.0);
    assert_ulps_eq!(grid_fn.y_max(), 10.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_y_range_single_extremum() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((-2.0, 2.0), 100, |x| -(x * x))?;
    assert!(grid_fn.y_max() > -0.001 && grid_fn.y_max() <= 0.0);
    assert!(grid_fn.y_min() < -3.9);
    Ok(())
  }

  #[test]
  fn test_gridfn_to_pairs() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0], array![10.0, 20.0, 30.0])?;
    let expected = vec![(0.0, 10.0), (1.0, 20.0), (2.0, 30.0)];
    let actual = grid_fn.to_pairs();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_to_pairs_roundtrip() -> Result<(), Report> {
    let pairs = vec![(0.0, 5.0), (1.0, 10.0), (2.0, 15.0), (3.0, 20.0)];
    let grid_fn = GridFn::from_pairs(pairs.clone())?;
    let actual = grid_fn.to_pairs();
    assert_eq!(pairs, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_to_pairs_nonuniform() -> Result<(), Report> {
    let pairs = vec![(0.0, 5.0), (0.1, 10.0), (1.0, 15.0), (10.0, 20.0)];
    let grid_fn = GridFn::from_pairs(pairs.clone())?;
    let actual = grid_fn.to_pairs();
    assert_eq!(pairs, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_clone() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0], array![0.0, 10.0])?;
    let cloned = grid_fn.clone();
    assert_eq!(grid_fn, cloned);
    Ok(())
  }

  #[test]
  fn test_gridfn_eq() -> Result<(), Report> {
    let grid_fn1 = GridFn::new(array![0.0, 1.0], array![0.0, 10.0])?;
    let grid_fn2 = GridFn::new(array![0.0, 1.0], array![0.0, 10.0])?;
    let grid_fn3 = GridFn::new(array![0.0, 1.0], array![0.0, 11.0])?;
    assert_eq!(grid_fn1, grid_fn2);
    assert_ne!(grid_fn1, grid_fn3);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many_large_uniform_spacing() -> Result<(), Report> {
    let n = 1000;
    let x = Array1::linspace(0.0, 100.0, n);
    let y = x.mapv(|v| v * v);
    let grid_fn = GridFn::new(x, y)?;
    let queries = array![25.0, 50.0, 75.0];
    let expected = array![625.0, 2500.0, 5625.0];
    let actual = grid_fn.interp_many(&queries)?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_large_nonuniform_spacing() -> Result<(), Report> {
    let n = 500;
    let x: Array1<f64> = (0..n).map(|i| (i as f64).powi(2) * 0.01).collect();
    let y = x.mapv(|v| v.sqrt());
    let grid_fn = GridFn::new(x, y)?;
    let query = 100.0;
    let result = grid_fn.interp(query)?;
    assert!(result > 9.9 && result < 10.1);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many_with_extrapolation() -> Result<(), Report> {
    let grid_fn = GridFn::new(Array1::linspace(0.0, 10.0, 100), Array1::linspace(0.0, 100.0, 100))?;
    let queries = Array1::linspace(-1.0, 11.0, 50);
    let expected = Array1::linspace(-10.0, 110.0, 50);
    let actual = grid_fn.interp_many(&queries)?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_negative_domain() -> Result<(), Report> {
    let grid_fn = GridFn::new(
      array![-10.0, -5.0, 0.0, 5.0, 10.0],
      array![100.0, 25.0, 0.0, 25.0, 100.0],
    )?;
    assert_ulps_eq!(grid_fn.interp(-7.5)?, 62.5);
    assert_ulps_eq!(grid_fn.interp(2.5)?, 12.5);
    assert_ulps_eq!(grid_fn.interp(-15.0)?, 175.0);
    assert_ulps_eq!(grid_fn.interp(15.0)?, 175.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_negative_values() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0, 3.0], array![10.0, -5.0, -20.0, 15.0])?;
    assert_ulps_eq!(grid_fn.interp(0.5)?, 2.5);
    assert_ulps_eq!(grid_fn.interp(1.5)?, -12.5);
    assert_ulps_eq!(grid_fn.interp(2.5)?, -2.5);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_oscillating_function() -> Result<(), Report> {
    let x = Array1::linspace(0.0, 10.0, 100);
    let y = x.mapv(|v: f64| (v * 2.0).sin());
    let grid_fn = GridFn::new(x, y)?;
    let result = grid_fn.interp(5.0)?;
    assert_abs_diff_eq!(result, (10.0_f64).sin(), epsilon = 0.01);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_steep_gradient() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 0.001, 1.0], array![0.0, 1000.0, 1001.0])?;
    assert_ulps_eq!(grid_fn.interp(0.0005)?, 500.0);
    assert_ulps_eq!(grid_fn.interp(0.5005)?, 1000.5, max_ulps = 1000000);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_boundary_precision() -> Result<(), Report> {
    let x = Array1::linspace(0.0, 1.0, 1000);
    let y = Array1::zeros(1000);
    let grid_fn = GridFn::new(x, y)?;
    assert_ulps_eq!(grid_fn.interp(0.0)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(1.0)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(0.999999)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(0.000001)?, 0.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many_at_boundaries() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0, 2.0, 3.0, 4.0], array![0.0, 1.0, 4.0, 9.0, 16.0])?;
    let queries = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 1.0, 4.0, 9.0, 16.0];
    let actual = grid_fn.interp_many(&queries)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_wide_range() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![1e-10, 1e10], array![0.0, 1.0])?;
    assert_ulps_eq!(grid_fn.interp(5e9)?, 0.5, max_ulps = 10000);
    assert_ulps_eq!(grid_fn.interp(1e-10)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(1e10)?, 1.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_piecewise_linear_segments() -> Result<(), Report> {
    let grid_fn = GridFn::new(
      array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
      array![0.0, 10.0, 10.0, 30.0, 30.0, 50.0],
    )?;
    assert_ulps_eq!(grid_fn.interp(0.5)?, 5.0);
    assert_ulps_eq!(grid_fn.interp(1.5)?, 10.0);
    assert_ulps_eq!(grid_fn.interp(2.5)?, 20.0);
    assert_ulps_eq!(grid_fn.interp(3.5)?, 30.0);
    assert_ulps_eq!(grid_fn.interp(4.5)?, 40.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_extrapolation_far_from_bounds() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![0.0, 1.0], array![0.0, 10.0])?;
    assert_ulps_eq!(grid_fn.interp(-100.0)?, -1000.0);
    assert_ulps_eq!(grid_fn.interp(100.0)?, 1000.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_zero_crossings() -> Result<(), Report> {
    let grid_fn = GridFn::new(array![-2.0, -1.0, 0.0, 1.0, 2.0], array![4.0, 1.0, 0.0, 1.0, 4.0])?;
    assert_ulps_eq!(grid_fn.interp(-0.5)?, 0.5);
    assert_ulps_eq!(grid_fn.interp(0.0)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(0.5)?, 0.5);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_grid_with_interp_exponential() -> Result<(), Report> {
    let grid_fn = GridFn::from_grid((0.0, 10.0), 0.01, |x| x.exp())?;
    assert_eq!(grid_fn.x().len(), 1001);
    assert_ulps_eq!(grid_fn.interp(0.0)?, 1.0);
    assert_ulps_eq!(grid_fn.interp(1.0)?, 1.0_f64.exp(), max_ulps = 100000);
    assert_ulps_eq!(grid_fn.interp(5.0)?, 5.0_f64.exp(), max_ulps = 10000000);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_n_points_quadratic() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 2.0), 5, |x| x * x)?;
    assert_eq!(grid_fn.x().len(), 5);
    assert_ulps_eq!(grid_fn.x()[0], 0.0);
    assert_ulps_eq!(grid_fn.x()[4], 2.0);
    assert_ulps_eq!(grid_fn.y()[0], 0.0);
    assert_ulps_eq!(grid_fn.y()[2], 1.0);
    assert_ulps_eq!(grid_fn.y()[4], 4.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_n_points_sine_with_interp() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, std::f64::consts::PI), 100, |x: f64| x.sin())?;
    assert_eq!(grid_fn.x().len(), 100);
    assert_ulps_eq!(grid_fn.interp(0.0)?, 0.0, epsilon = 1e-10);
    assert_ulps_eq!(grid_fn.interp(std::f64::consts::PI / 2.0)?, 1.0, epsilon = 1e-2);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_fn_samples_uniform() -> Result<(), Report> {
    let x = Array1::linspace(0.0, 1.0, 11);
    let grid_fn = GridFn::from_fn_samples(x, |x| x * x * x)?;
    assert_eq!(grid_fn.x().len(), 11);
    assert_ulps_eq!(grid_fn.y()[0], 0.0);
    assert_ulps_eq!(grid_fn.y()[5], 0.125);
    assert_ulps_eq!(grid_fn.y()[10], 1.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_fn_samples_nonuniform() -> Result<(), Report> {
    let x = array![0.0, 0.1, 0.5, 1.0, 2.0, 5.0];
    let grid_fn = GridFn::from_fn_samples(x, |x| x.ln())?;
    assert_eq!(grid_fn.x().len(), 6);
    assert_ulps_eq!(grid_fn.y()[3], 0.0_f64);
    assert_ulps_eq!(grid_fn.y()[5], 5.0_f64.ln());
    Ok(())
  }

  #[test]
  fn test_gridfn_constant() -> Result<(), Report> {
    let grid_fn = GridFn::constant((0.0, 10.0), 20, 42.0)?;
    assert_eq!(grid_fn.x().len(), 20);
    assert_eq!(grid_fn.y().len(), 20);
    assert_ulps_eq!(grid_fn.y()[0], 42.0);
    assert_ulps_eq!(grid_fn.y()[10], 42.0);
    assert_ulps_eq!(grid_fn.y()[19], 42.0);
    assert_ulps_eq!(grid_fn.interp(5.0)?, 42.0);
    assert_ulps_eq!(grid_fn.interp(-5.0)?, 42.0);
    assert_ulps_eq!(grid_fn.interp(15.0)?, 42.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_zeros() -> Result<(), Report> {
    let grid_fn = GridFn::zeros((-1.0, 1.0), 10)?;
    assert_eq!(grid_fn.x().len(), 10);
    assert_eq!(grid_fn.y().len(), 10);
    assert_ulps_eq!(grid_fn.y()[0], 0.0);
    assert_ulps_eq!(grid_fn.y()[5], 0.0);
    assert_ulps_eq!(grid_fn.y()[9], 0.0);
    assert_ulps_eq!(grid_fn.interp(0.0)?, 0.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_ones() -> Result<(), Report> {
    let grid_fn = GridFn::ones((0.0, 5.0), 15)?;
    assert_eq!(grid_fn.x().len(), 15);
    assert_eq!(grid_fn.y().len(), 15);
    assert_ulps_eq!(grid_fn.y()[0], 1.0);
    assert_ulps_eq!(grid_fn.y()[7], 1.0);
    assert_ulps_eq!(grid_fn.y()[14], 1.0);
    assert_ulps_eq!(grid_fn.interp(2.5)?, 1.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_pairs() -> Result<(), Report> {
    let pairs = vec![(0.0, 0.0), (1.0, 10.0), (2.0, 20.0), (3.0, 30.0)];
    let grid_fn = GridFn::from_pairs(pairs)?;
    assert_eq!(grid_fn.x().len(), 4);
    assert_eq!(grid_fn.y().len(), 4);
    assert_ulps_eq!(grid_fn.x()[0], 0.0);
    assert_ulps_eq!(grid_fn.x()[3], 3.0);
    assert_ulps_eq!(grid_fn.y()[0], 0.0);
    assert_ulps_eq!(grid_fn.y()[3], 30.0);
    assert_ulps_eq!(grid_fn.interp(1.5)?, 15.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_from_pairs_nonuniform() -> Result<(), Report> {
    let pairs = vec![(0.0, 5.0), (0.1, 10.0), (1.0, 15.0), (10.0, 20.0)];
    let grid_fn = GridFn::from_pairs(pairs)?;
    assert_eq!(grid_fn.x().len(), 4);
    assert_ulps_eq!(grid_fn.interp(0.05)?, 7.5);
    assert_ulps_eq!(grid_fn.interp(5.5)?, 17.5);
    Ok(())
  }
}
