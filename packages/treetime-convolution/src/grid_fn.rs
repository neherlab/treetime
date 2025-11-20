use std::cmp::min;
use std::fmt::Debug;

use crate::InterpElem;
use crate::grid::Grid;
use crate::grid_iter::GridIter;
use crate::interp_nonuniform::interp_nonuniform;
use approx::UlpsEq;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, s};
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_utils::make_error;
use treetime_utils::ndarray::has_uniform_spacing;
use treetime_utils::serde::{array1_as_vec, array1_from_vec};

impl<T: InterpElem> IntoIterator for &Grid<T>
where
  T: Float,
{
  type Item = T;
  type IntoIter = GridIter<T>;

  fn into_iter(self) -> Self::IntoIter {
    self.iter()
  }
}

/// Function represented on a uniform grid for piecewise linear interpolation
///
/// Represents a function as a set of (x, y) points on a uniformly-spaced grid, providing
/// linear interpolation between points and linear extrapolation beyond the grid boundaries.
///
/// # Invariants
///
/// - Grid must be uniformly spaced with spacing `dx`
/// - Grid starts at `x_min`
/// - `y` array must contain at least 2 points
/// - `dx` must be positive
/// - These invariants are enforced by the type system and cannot be violated
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound(serialize = "T: Serialize", deserialize = "T: Deserialize<'de>"))]
pub struct GridFn<T: InterpElem> {
  grid: Grid<T>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  y: Array1<T>,
}

impl<T: InterpElem> GridFn<T> {
  pub fn from_grid_array(grid: Grid<T>, y: Array1<T>) -> Result<Self, Report> {
    if grid.n_points() != y.len() {
      return make_error!(
        "Grid has {} points but y array has {} elements",
        grid.n_points(),
        y.len()
      );
    }
    Ok(Self { grid, y })
  }

  pub fn from_grid_fn<F>(grid: Grid<T>, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let n_points = grid.n_points();
    let y = Array1::from_shape_fn(n_points, |i| y_fn(grid.x_at(i)));
    Self::from_grid_array(grid, y)
  }

  pub fn from_arrays(x: &Array1<T>, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    if x.len() != y.len() {
      return make_error!(
        "x and y arrays must have the same length, got {} and {}",
        x.len(),
        y.len()
      );
    }
    let grid = Grid::from_array(x)?;
    Self::from_grid_array(grid, y)
  }

  /// Constructs GridFn from non-uniformly spaced arrays by resampling to uniform grid
  ///
  /// Takes non-uniform (x, y) arrays and resamples them to a uniform grid using linear
  /// interpolation. The grid spacing is determined by the smallest spacing in the input
  /// to preserve detail.
  ///
  /// # Arguments
  ///
  /// * `x` - Non-uniform x coordinates (must be sorted ascending)
  /// * `y` - Corresponding y values
  ///
  /// # Returns
  ///
  /// GridFn with uniformly spaced grid covering the same range as input
  pub fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    if x.len() < 2 {
      return make_error!("Grid must have at least 2 points, got {}", x.len());
    }
    if x.len() != y.len() {
      return make_error!(
        "x and y arrays must have the same length, got {} and {}",
        x.len(),
        y.len()
      );
    }

    if has_uniform_spacing(x) {
      let grid = Grid::from_array(x)?;
      return Self::from_grid_array(grid, y.clone());
    }

    let x_min = x[0];
    let x_max = x[x.len() - 1];
    let dx = find_min_spacing(x)?;
    let n_points = ((x_max - x_min) / dx).ceil().to_usize().unwrap() + 1;
    let dx = (x_max - x_min) / T::from(n_points - 1).unwrap();

    if n_points > 1_000_000 {
      return make_error!("Resampling would require {n_points} points, which exceeds safety limit");
    }

    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    let y_uniform = interp_nonuniform(x, y, n_points, |i| grid.x_at(i))?;
    Self::from_grid_array(grid, y_uniform)
  }

  pub fn from_n_points<F>((x_min, x_max): (T, T), n_points: usize, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    Self::from_grid_fn(grid, y_fn)
  }

  pub fn from_grid<F>((x_min, x_max): (T, T), dx: T, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid = Grid::from_range_dx(x_min, x_max, dx)?;
    Self::from_grid_fn(grid, y_fn)
  }

  pub fn from_start_dx_values(x_min: T, dx: T, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_start_dx(x_min, dx, y.len())?;
    Self::from_grid_array(grid, y)
  }

  pub fn from_range_values((x_min, x_max): (T, T), y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, y.len())?;
    Self::from_grid_array(grid, y)
  }

  pub fn constant((x_min, x_max): (T, T), n_points: usize, value: T) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    let y = Array1::from_elem(n_points, value);
    Self::from_grid_array(grid, y)
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

  // TODO: inefficient. Try to remove this method
  pub fn x(&self) -> Array1<T>
  where
    T: Float,
  {
    self.grid.to_array()
  }

  pub fn y(&self) -> &Array1<T> {
    &self.y
  }

  pub fn x_min(&self) -> T {
    self.grid.x_min()
  }

  pub fn x_max(&self) -> T
  where
    T: Float,
  {
    self.grid.x_max()
  }

  pub fn x_range(&self) -> (T, T)
  where
    T: Float,
  {
    self.grid.x_range()
  }

  pub fn n_points(&self) -> usize {
    self.grid.n_points()
  }

  pub fn len(&self) -> usize {
    self.grid.len()
  }

  pub fn is_empty(&self) -> bool {
    self.grid.is_empty()
  }

  pub fn dx(&self) -> T {
    self.grid.dx()
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

  pub fn to_pairs(&self) -> Vec<(T, T)>
  where
    T: Float,
  {
    izip!(self.grid.iter(), self.y.iter())
      .map(|(x, &y)| (x, y))
      .collect_vec()
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
  pub fn interp(&self, xi: T) -> Result<T, Report>
  where
    T: Float,
  {
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();

    if xi < x_min {
      return Ok(self.extrapolate_left(xi));
    }

    if xi > x_max {
      return Ok(self.extrapolate_right(xi));
    }

    let idx = self.grid.find_interval_index(xi);
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
  pub fn interp_many(&self, queries: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float,
  {
    debug_assert!(
      queries.windows(2).into_iter().all(|w| w[0] <= w[1]),
      "queries must be sorted in non-decreasing order"
    );

    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();
    let n = self.grid.n_points();
    let mut result = Array1::zeros(queries.len());
    let mut grid_idx = 0_usize;

    for (i, &query) in queries.iter().enumerate() {
      if query < x_min {
        result[i] = self.extrapolate_left(query);
      } else if query > x_max {
        result[i] = self.extrapolate_right(query);
      } else {
        while grid_idx < n - 1 && self.grid.x_at(grid_idx + 1) < query {
          grid_idx += 1;
        }
        let idx = min(grid_idx, n - 2);
        result[i] = self.interpolate_at(query, idx);
      }
    }
    Ok(result)
  }

  fn extrapolate_left(&self, q: T) -> T
  where
    T: Float,
  {
    let slope = (self.y[1] - self.y[0]) / self.grid.dx();
    self.y[0] + slope * (q - self.grid.x_min())
  }

  fn extrapolate_right(&self, q: T) -> T
  where
    T: Float,
  {
    let n = self.grid.n_points();
    let slope = (self.y[n - 1] - self.y[n - 2]) / self.grid.dx();
    let x_max = self.grid.x_max();
    self.y[n - 1] + slope * (q - x_max)
  }

  fn interpolate_at(&self, q: T, idx: usize) -> T
  where
    T: Float,
  {
    let n = self.grid.n_points();
    if idx >= n - 1 {
      return self.y[n - 1];
    }
    let x0 = self.grid.x_at(idx);
    let y0 = self.y[idx];
    let y1 = self.y[idx + 1];
    let t = (q - x0) / self.grid.dx();
    y0 + t * (y1 - y0)
  }

  #[must_use]
  pub fn mapv<F>(&self, f: F) -> Self
  where
    F: Fn(T) -> T,
  {
    Self {
      grid: self.grid,
      y: self.y.mapv(f),
    }
  }

  pub fn mapv_inplace<F>(&mut self, f: F)
  where
    F: Fn(T) -> T,
  {
    self.y.mapv_inplace(f);
  }

  /// Negates the argument of the function: f(x) -> f(-x).
  /// This reflects the function across the y-axis.
  /// The domain [x_min, x_max] becomes [-x_max, -x_min].
  /// The y-values are reversed.
  #[must_use]
  pub fn negate_arg(&self) -> Self
  where
    T: Float,
  {
    let mut result = self.clone();
    result.negate_arg_inplace();
    result
  }

  /// Negates the argument of the function in-place: f(x) -> f(-x).
  /// This reflects the function across the y-axis.
  /// The domain [x_min, x_max] becomes [-x_max, -x_min].
  /// The y-values are reversed.
  pub fn negate_arg_inplace(&mut self)
  where
    T: Float,
  {
    let x_max = self.grid.x_max();
    let dx = self.grid.dx();
    let n_points = self.grid.n_points();
    self.grid = Grid::from_start_dx(-x_max, dx, n_points).unwrap();

    // Reverse y in-place
    let n = self.y.len();
    for i in 0..n / 2 {
      self.y.swap(i, n - 1 - i);
    }
  }

  /// Resamples function to a new uniform grid with specified spacing
  ///
  /// Creates a new GridFn on a uniform grid with the given spacing, using linear
  /// interpolation (and extrapolation if needed) from the current function.
  ///
  /// # Arguments
  ///
  /// * `x_range` - New grid range (x_min, x_max)
  /// * `dx` - Grid spacing for the new grid
  ///
  /// # Returns
  ///
  /// New GridFn with uniformly spaced grid covering the specified range
  pub fn resample_to_grid(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_dx(x_min, x_max, dx)?;
    let n_points = grid.n_points();
    let mut y_new = Array1::zeros(n_points);
    for i in 0..n_points {
      let x = grid.x_at(i);
      y_new[i] = self.interp(x)?;
    }
    Self::from_grid_array(grid, y_new)
  }

  /// Resamples function to a new uniform grid with specified number of points
  ///
  /// Creates a new GridFn on a uniform grid with the given number of points, using
  /// linear interpolation (and extrapolation if needed) from the current function.
  ///
  /// # Arguments
  ///
  /// * `x_range` - New grid range (x_min, x_max)
  /// * `n_points` - Number of points in the new grid
  ///
  /// # Returns
  ///
  /// New GridFn with uniformly spaced grid covering the specified range
  pub fn resample_to_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    let mut y_new = Array1::zeros(n_points);
    for i in 0..n_points {
      let x = grid.x_at(i);
      y_new[i] = self.interp(x)?;
    }
    Self::from_grid_array(grid, y_new)
  }
}

/// Finds the smallest spacing between consecutive points in a sorted array
///
/// Used to determine optimal grid spacing when resampling non-uniform data.
/// Preserves maximum detail by using the finest resolution present in the input.
///
/// # Arguments
///
/// * `x` - Sorted array of x coordinates (must be ascending)
///
/// # Returns
///
/// Minimum spacing between consecutive points
fn find_min_spacing<T>(x: &Array1<T>) -> Result<T, Report>
where
  T: Float + Debug,
{
  if x.len() < 2 {
    return make_error!("Array must have at least 2 points");
  }

  let diffs = &x.slice(s![1..]) - &x.slice(s![..-1]);

  if diffs.iter().any(|&dx| dx <= T::zero()) {
    return make_error!("x array must be sorted in ascending order");
  }

  let min_dx = *diffs.min().unwrap();

  if !min_dx.is_finite() || min_dx <= T::zero() {
    return make_error!("Invalid spacing in input array: {min_dx:?}");
  }

  Ok(min_dx)
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
  #[case((0.0, 1.0), array![0.0, 10.0], 0.0, 0.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], 1.0, 10.0)]
  #[case((0.0, 2.0), array![0.0, 5.0, 10.0], 1.0, 5.0)]
  #[case((-5.0, 5.0), array![100.0, 200.0], -5.0, 100.0)]
  #[case((-5.0, 5.0), array![100.0, 200.0], 5.0, 200.0)]
  #[trace]
  fn test_gridfn_interp_exact_grid_points(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case((0.0, 1.0), array![0.0, 10.0], 0.5, 5.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], 0.3, 3.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], 0.7, 7.0)]
  #[case((0.0, 2.0), array![10.0, 20.0], 1.0, 15.0)]
  #[case((-1.0, 1.0), array![0.0, 100.0], 0.0, 50.0)]
  #[case((0.0, 2.0), array![0.0, 10.0, 20.0], 0.5, 5.0)]
  #[case((0.0, 2.0), array![0.0, 10.0, 20.0], 1.5, 15.0)]
  #[trace]
  fn test_gridfn_interp_interior_points(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case((0.0, 1.0), array![0.0, 10.0], -1.0, -10.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], -0.5, -5.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], -2.0, -20.0)]
  #[case((1.0, 2.0), array![5.0, 15.0], 0.0, -5.0)]
  #[case((1.0, 3.0), array![10.0, 20.0, 30.0], 0.0, 0.0)]
  #[trace]
  fn test_gridfn_interp_left_extrapolation(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case((0.0, 1.0), array![0.0, 10.0], 2.0, 20.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], 1.5, 15.0)]
  #[case((0.0, 1.0), array![0.0, 10.0], 3.0, 30.0)]
  #[case((1.0, 2.0), array![5.0, 15.0], 3.0, 25.0)]
  #[case((0.0, 2.0), array![10.0, 20.0, 30.0], 3.0, 40.0)]
  #[trace]
  fn test_gridfn_interp_right_extrapolation(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case((0.0, 1.0), array![5.0, 5.0], 0.5, 5.0)]
  #[case((0.0, 1.0), array![5.0, 5.0], -1.0, 5.0)]
  #[case((0.0, 1.0), array![5.0, 5.0], 2.0, 5.0)]
  #[case((0.0, 2.0), array![42.0, 42.0, 42.0], 1.5, 42.0)]
  #[trace]
  fn test_gridfn_interp_constant_function(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-14);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case((0.0, 1.0), array![0.0, 1e-15], 0.5, 5e-16)]
  #[case((0.0, 1e-10), array![0.0, 1.0], 5e-11, 0.5)]
  #[case((1e10, 1e10 + 1.0), array![0.0, 1.0], 1e10 + 0.5, 0.5)]
  #[trace]
  fn test_gridfn_interp_numeric_precision(
    #[case] x_range: (f64, f64),
    #[case] y: Array1<f64>,
    #[case] query: f64,
    #[case] expected: f64,
  ) -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    let actual = grid_fn.interp(query)?;
    assert_ulps_eq!(expected, actual, epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
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
    let x_range = (0.0, 2.0);
    let y = array![10.0, 20.0, 30.0];
    let grid_fn = GridFn::from_range_values(x_range, y.clone())?;
    let x = grid_fn.x();
    assert_ulps_eq!(x[0], 0.0);
    assert_ulps_eq!(x[1], 1.0);
    assert_ulps_eq!(x[2], 2.0);
    assert_eq!(&y, grid_fn.y());
    Ok(())
  }

  #[test]
  fn test_gridfn_x_min_x_max() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 3.0), array![10.0, 20.0, 30.0, 40.0])?;
    assert_ulps_eq!(grid_fn.x_min(), 0.0);
    assert_ulps_eq!(grid_fn.x_max(), 3.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_x_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((-5.0, 5.0), array![100.0, 50.0, 0.0])?;
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
    let grid_fn = GridFn::from_range_values((0.0, 3.0), array![10.0, 5.0, 30.0, 15.0])?;
    assert_ulps_eq!(grid_fn.y_min(), 5.0);
    assert_ulps_eq!(grid_fn.y_max(), 30.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_y_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 3.0), array![-10.0, 20.0, 5.0, -5.0])?;
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
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![10.0, 20.0, 30.0])?;
    let expected = vec![(0.0, 10.0), (1.0, 20.0), (2.0, 30.0)];
    let actual = grid_fn.to_pairs();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_to_pairs_roundtrip() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 3.0), array![5.0, 10.0, 15.0, 20.0])?;
    let expected = vec![(0.0, 5.0), (1.0, 10.0), (2.0, 15.0), (3.0, 20.0)];
    let actual = grid_fn.to_pairs();
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_clone() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let cloned = grid_fn.clone();
    assert_eq!(grid_fn, cloned);
    Ok(())
  }

  #[test]
  fn test_gridfn_eq() -> Result<(), Report> {
    let grid_fn1 = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let grid_fn2 = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let grid_fn3 = GridFn::from_range_values((0.0, 1.0), array![0.0, 11.0])?;
    assert_eq!(grid_fn1, grid_fn2);
    assert_ne!(grid_fn1, grid_fn3);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many_large_uniform_spacing() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 100.0), 1000, |v| v * v)?;
    let queries = array![25.0, 50.0, 75.0];
    let expected = array![625.0, 2500.0, 5625.0];
    let actual = grid_fn.interp_many(&queries)?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many_with_extrapolation() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 10.0), 100, |v| v * 10.0)?;
    let queries = Array1::linspace(-1.0, 11.0, 50);
    let expected = Array1::linspace(-10.0, 110.0, 50);
    let actual = grid_fn.interp_many(&queries)?;
    assert_abs_diff_eq!(expected, actual, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_negative_domain() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((-10.0, 10.0), array![100.0, 25.0, 0.0, 25.0, 100.0])?;
    assert_ulps_eq!(grid_fn.interp(-7.5)?, 62.5);
    assert_ulps_eq!(grid_fn.interp(2.5)?, 12.5);
    assert_ulps_eq!(grid_fn.interp(-15.0)?, 175.0);
    assert_ulps_eq!(grid_fn.interp(15.0)?, 175.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_negative_values() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 3.0), array![10.0, -5.0, -20.0, 15.0])?;
    assert_ulps_eq!(grid_fn.interp(0.5)?, 2.5);
    assert_ulps_eq!(grid_fn.interp(1.5)?, -12.5);
    assert_ulps_eq!(grid_fn.interp(2.5)?, -2.5);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_oscillating_function() -> Result<(), Report> {
    let grid_fn = GridFn::from_n_points((0.0, 10.0), 100, |v: f64| (v * 2.0).sin())?;
    let result = grid_fn.interp(5.0)?;
    assert_abs_diff_eq!(result, (10.0_f64).sin(), epsilon = 0.01);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_steep_gradient() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 1000.0, 1001.0])?;
    assert_ulps_eq!(grid_fn.interp(0.25)?, 500.0);
    assert_ulps_eq!(grid_fn.interp(0.75)?, 1000.5, max_ulps = 1000000);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_boundary_precision() -> Result<(), Report> {
    let grid_fn = GridFn::zeros((0.0, 1.0), 1000)?;
    assert_ulps_eq!(grid_fn.interp(0.0)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(1.0)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(0.999999)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(0.000001)?, 0.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_many_at_boundaries() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 4.0), array![0.0, 1.0, 4.0, 9.0, 16.0])?;
    let queries = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let expected = array![0.0, 1.0, 4.0, 9.0, 16.0];
    let actual = grid_fn.interp_many(&queries)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_wide_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((1e-10, 1e10), array![0.0, 1.0])?;
    assert_ulps_eq!(grid_fn.interp(5e9)?, 0.5, max_ulps = 10000);
    assert_ulps_eq!(grid_fn.interp(1e-10)?, 0.0);
    assert_ulps_eq!(grid_fn.interp(1e10)?, 1.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_piecewise_linear_segments() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 5.0), array![0.0, 10.0, 10.0, 30.0, 30.0, 50.0])?;
    assert_ulps_eq!(grid_fn.interp(0.5)?, 5.0);
    assert_ulps_eq!(grid_fn.interp(1.5)?, 10.0);
    assert_ulps_eq!(grid_fn.interp(2.5)?, 20.0);
    assert_ulps_eq!(grid_fn.interp(3.5)?, 30.0);
    assert_ulps_eq!(grid_fn.interp(4.5)?, 40.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_extrapolation_far_from_bounds() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    assert_ulps_eq!(grid_fn.interp(-100.0)?, -1000.0);
    assert_ulps_eq!(grid_fn.interp(100.0)?, 1000.0);
    Ok(())
  }

  #[test]
  fn test_gridfn_interp_zero_crossings() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((-2.0, 2.0), array![4.0, 1.0, 0.0, 1.0, 4.0])?;
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
  fn test_gridfn_mapv() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let expected = GridFn::from_range_values((0.0, 2.0), array![2.0, 4.0, 6.0])?;
    let actual = grid_fn.mapv(|y| y * 2.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_square() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 3.0), array![1.0, 2.0, 3.0, 4.0])?;
    let expected = GridFn::from_range_values((0.0, 3.0), array![1.0, 4.0, 9.0, 16.0])?;
    let actual = grid_fn.mapv(|y| y * y);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_negate() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, -2.0, 3.0])?;
    let expected = GridFn::from_range_values((0.0, 2.0), array![-1.0, 2.0, -3.0])?;
    let actual = grid_fn.mapv(|y| -y);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_preserves_grid() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((5.0, 10.0), array![1.0, 2.0, 3.0])?;
    let mapped = grid_fn.mapv(|y| y + 10.0);
    assert_ulps_eq!(mapped.x_min(), grid_fn.x_min());
    assert_ulps_eq!(mapped.x_max(), grid_fn.x_max());
    assert_ulps_eq!(mapped.dx(), grid_fn.dx());
    assert_eq!(mapped.n_points(), grid_fn.n_points());
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_chain() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let expected = GridFn::from_range_values((0.0, 2.0), array![4.0, 6.0, 8.0])?;
    let actual = grid_fn.mapv(|y| y * 2.0).mapv(|y| y + 2.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_inplace() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    grid_fn.mapv_inplace(|y| y * 2.0);
    let expected = GridFn::from_range_values((0.0, 2.0), array![2.0, 4.0, 6.0])?;
    assert_eq!(expected, grid_fn);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_inplace_square() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((0.0, 3.0), array![1.0, 2.0, 3.0, 4.0])?;
    grid_fn.mapv_inplace(|y| y * y);
    let expected = GridFn::from_range_values((0.0, 3.0), array![1.0, 4.0, 9.0, 16.0])?;
    assert_eq!(expected, grid_fn);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_inplace_negate() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, -2.0, 3.0])?;
    grid_fn.mapv_inplace(|y| -y);
    let expected = GridFn::from_range_values((0.0, 2.0), array![-1.0, 2.0, -3.0])?;
    assert_eq!(expected, grid_fn);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_inplace_preserves_grid() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((5.0, 10.0), array![1.0, 2.0, 3.0])?;
    let x_min = grid_fn.x_min();
    let x_max = grid_fn.x_max();
    let dx = grid_fn.dx();
    let n_points = grid_fn.n_points();
    grid_fn.mapv_inplace(|y| y + 10.0);
    assert_ulps_eq!(grid_fn.x_min(), x_min);
    assert_ulps_eq!(grid_fn.x_max(), x_max);
    assert_ulps_eq!(grid_fn.dx(), dx);
    assert_eq!(grid_fn.n_points(), n_points);
    Ok(())
  }

  #[test]
  fn test_gridfn_mapv_inplace_chain() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    grid_fn.mapv_inplace(|y| y * 2.0);
    grid_fn.mapv_inplace(|y| y + 2.0);
    let expected = GridFn::from_range_values((0.0, 2.0), array![4.0, 6.0, 8.0])?;
    assert_eq!(expected, grid_fn);
    Ok(())
  }

  #[test]
  fn test_gridfn_negate_arg_inplace() -> Result<(), Report> {
    let mut grid_fn = GridFn::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    grid_fn.negate_arg_inplace();
    let expected = GridFn::from_range_values((-2.0, 0.0), array![3.0, 2.0, 1.0])?;
    assert_eq!(expected, grid_fn);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_grid_same_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 4.0), array![0.0, 1.0, 4.0, 9.0, 16.0])?;
    let resampled = grid_fn.resample_to_grid((0.0, 4.0), 1.0)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 4.0);
    assert_ulps_eq!(resampled.dx(), 1.0);
    assert_eq!(resampled.n_points(), 5);
    assert_abs_diff_eq!(resampled.y(), &array![0.0, 1.0, 4.0, 9.0, 16.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_grid_finer() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
    let resampled = grid_fn.resample_to_grid((0.0, 2.0), 0.5)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 2.0);
    assert_ulps_eq!(resampled.dx(), 0.5);
    assert_eq!(resampled.n_points(), 5);
    assert_abs_diff_eq!(resampled.y(), &array![0.0, 5.0, 10.0, 15.0, 20.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_grid_coarser() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 4.0), array![0.0, 1.0, 4.0, 9.0, 16.0])?;
    let resampled = grid_fn.resample_to_grid((0.0, 4.0), 2.0)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 4.0);
    assert_ulps_eq!(resampled.dx(), 2.0);
    assert_eq!(resampled.n_points(), 3);
    assert_abs_diff_eq!(resampled.y(), &array![0.0, 4.0, 16.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_grid_different_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
    let resampled = grid_fn.resample_to_grid((1.0, 3.0), 0.5)?;
    assert_ulps_eq!(resampled.x_min(), 1.0);
    assert_ulps_eq!(resampled.x_max(), 3.0);
    assert_ulps_eq!(resampled.dx(), 0.5);
    assert_eq!(resampled.n_points(), 5);
    assert_abs_diff_eq!(resampled.y(), &array![10.0, 15.0, 20.0, 25.0, 30.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_grid_with_extrapolation() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let resampled = grid_fn.resample_to_grid((-1.0, 2.0), 0.5)?;
    assert_ulps_eq!(resampled.x_min(), -1.0);
    assert_ulps_eq!(resampled.x_max(), 2.0);
    assert_eq!(resampled.n_points(), 7);
    assert_abs_diff_eq!(
      resampled.y(),
      &array![-10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0],
      epsilon = 1e-10
    );
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_n_points_same_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 4.0), array![0.0, 1.0, 4.0, 9.0, 16.0])?;
    let resampled = grid_fn.resample_to_n_points((0.0, 4.0), 5)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 4.0);
    assert_eq!(resampled.n_points(), 5);
    assert_abs_diff_eq!(resampled.y(), &array![0.0, 1.0, 4.0, 9.0, 16.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_n_points_finer() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
    let resampled = grid_fn.resample_to_n_points((0.0, 2.0), 5)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 2.0);
    assert_eq!(resampled.n_points(), 5);
    assert_abs_diff_eq!(resampled.y(), &array![0.0, 5.0, 10.0, 15.0, 20.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_n_points_coarser() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 4.0), array![0.0, 1.0, 4.0, 9.0, 16.0])?;
    let resampled = grid_fn.resample_to_n_points((0.0, 4.0), 3)?;
    assert_ulps_eq!(resampled.x_min(), 0.0);
    assert_ulps_eq!(resampled.x_max(), 4.0);
    assert_eq!(resampled.n_points(), 3);
    assert_abs_diff_eq!(resampled.y(), &array![0.0, 4.0, 16.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_n_points_different_range() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 2.0), array![0.0, 10.0, 20.0])?;
    let resampled = grid_fn.resample_to_n_points((1.0, 3.0), 5)?;
    assert_ulps_eq!(resampled.x_min(), 1.0);
    assert_ulps_eq!(resampled.x_max(), 3.0);
    assert_eq!(resampled.n_points(), 5);
    assert_abs_diff_eq!(resampled.y(), &array![10.0, 15.0, 20.0, 25.0, 30.0], epsilon = 1e-10);
    Ok(())
  }

  #[test]
  fn test_gridfn_resample_to_n_points_with_extrapolation() -> Result<(), Report> {
    let grid_fn = GridFn::from_range_values((0.0, 1.0), array![0.0, 10.0])?;
    let resampled = grid_fn.resample_to_n_points((-1.0, 2.0), 7)?;
    assert_ulps_eq!(resampled.x_min(), -1.0);
    assert_ulps_eq!(resampled.x_max(), 2.0);
    assert_eq!(resampled.n_points(), 7);
    assert_abs_diff_eq!(
      resampled.y(),
      &array![-10.0, -5.0, 0.0, 5.0, 10.0, 15.0, 20.0],
      epsilon = 1e-10
    );
    Ok(())
  }
}
