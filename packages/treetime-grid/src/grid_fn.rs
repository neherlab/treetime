use std::fmt::Debug;

use crate::InterpElem;
use crate::grid::Grid;
use crate::interp_nonuniform::interp_nonuniform;
use approx::UlpsEq;
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::{Array1, s};
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_utils::array::ndarray::has_uniform_spacing;
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};
use treetime_utils::make_error;

/// Function represented on a uniform grid for piecewise linear interpolation
///
/// Represents a function as a set of (x, y) points on a uniformly-spaced grid, providing
/// linear interpolation between points and constant extrapolation beyond the grid boundaries.
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

  pub fn grid(&self) -> &Grid<T> {
    &self.grid
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
  /// Uses piecewise linear interpolation within grid bounds and constant extrapolation
  /// outside bounds (returns boundary values).
  ///
  /// # Arguments
  ///
  /// * `xi` - Point at which to evaluate the function
  ///
  /// # Returns
  ///
  /// Interpolated or extrapolated function value at `xi`
  pub fn interp(&self, xi: T) -> T
  where
    T: Float,
  {
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();

    if xi < x_min {
      return self.extrapolate_left(xi);
    }

    if xi > x_max {
      return self.extrapolate_right(xi);
    }

    let idx = self.grid.find_interval_index(xi);
    self.interpolate_at(xi, idx)
  }

  /// Interpolate function values at multiple points
  ///
  /// Uses piecewise linear interpolation within grid bounds and constant extrapolation
  /// outside bounds.
  ///
  /// # Arguments
  ///
  /// * `queries` - Query points at which to evaluate the function
  ///
  /// # Returns
  ///
  /// Array of interpolated or extrapolated function values at each query point
  pub fn interp_many(&self, queries: &Array1<T>) -> Array1<T>
  where
    T: Float,
  {
    queries.mapv(|q| self.interp(q))
  }

  fn extrapolate_left(&self, _q: T) -> T
  where
    T: Float,
  {
    // Use constant extrapolation: return first value
    // This matches Python scipy.interpolate.interp1d behavior when
    // the domain is extended with boundary values
    self.y[0]
  }

  fn extrapolate_right(&self, _q: T) -> T
  where
    T: Float,
  {
    // Use constant extrapolation: return last value
    // This matches Python scipy.interpolate.interp1d behavior when
    // the domain is extended with boundary values
    let n = self.grid.n_points();
    self.y[n - 1]
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

  /// Resamples function to a new grid
  ///
  /// Creates a new GridFn on the specified grid, using linear interpolation
  /// (and extrapolation if needed) from the current function.
  ///
  /// # Arguments
  ///
  /// * `grid` - Target grid for resampling
  ///
  /// # Returns
  ///
  /// New GridFn with values interpolated onto the target grid
  pub fn resample(&self, grid: &Grid<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let n_points = grid.n_points();
    let y_new = Array1::from_shape_fn(n_points, |i| self.interp(grid.x_at(i)));
    Self::from_grid_array(*grid, y_new)
  }

  /// Resamples function to a new uniform grid with specified start, spacing, and length
  ///
  /// Creates a new GridFn on a uniform grid defined by starting point, spacing,
  /// and number of points, using linear interpolation (and extrapolation if needed)
  /// from the current function.
  ///
  /// # Arguments
  ///
  /// * `x_min` - Starting x coordinate
  /// * `dx` - Grid spacing
  /// * `n_points` - Number of points in the new grid
  ///
  /// # Returns
  ///
  /// New GridFn with uniformly spaced grid
  pub fn resample_start_dx(&self, x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    self.resample(&grid)
  }

  /// Resamples function to a new uniform grid with specified range and number of points
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
  pub fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    self.resample(&grid)
  }

  /// Resamples function to a new uniform grid with specified range and spacing
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
  pub fn resample_range_dx(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_dx(x_min, x_max, dx)?;
    self.resample(&grid)
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
