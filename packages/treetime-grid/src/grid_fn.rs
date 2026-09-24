use std::fmt::Debug;

use crate::InterpElem;
use crate::boundary_behavior::BoundaryBehavior;
use crate::grid::Grid;
use crate::grid_edge::GridEdge;
use crate::hard_approach_law::Side;
use crate::interp_nonuniform::interp_nonuniform;
use approx::{UlpsEq, ulps_eq};
use eyre::Report;
use ndarray::{Array1, s};
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_utils::array::ndarray::has_uniform_spacing;
use treetime_utils::array::serde::{array1_as_vec, array1_from_vec};
use treetime_utils::make_error;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(bound(serialize = "T: Serialize", deserialize = "T: Deserialize<'de>"))]
pub struct GridFn<T: InterpElem> {
  grid: Grid<T>,
  #[serde(serialize_with = "array1_as_vec", deserialize_with = "array1_from_vec")]
  y: Array1<T>,
}

impl<T: InterpElem> GridFn<T> {
  fn from_grid_array(grid: Grid<T>, y: Array1<T>) -> Result<Self, Report> {
    if grid.n_points() != y.len() {
      return make_error!(
        "Grid has {} points but y array has {} elements",
        grid.n_points(),
        y.len()
      );
    }
    Ok(Self { grid, y })
  }

  fn from_grid_fn<F>(grid: Grid<T>, y_fn: F) -> Result<Self, Report>
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

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
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

  pub(crate) fn n_points(&self) -> usize {
    self.grid.n_points()
  }

  pub fn len(&self) -> usize {
    self.grid.len()
  }

  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  pub fn dx(&self) -> T {
    self.grid.dx()
  }

  pub fn y_min(&self) -> T {
    self
      .y
      .iter()
      .copied()
      .skip(1)
      .fold(self.y[0], |a, b| if a < b { a } else { b })
  }

  pub fn interp(&self, xi: T) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    self.interp_with_extrap(xi, BoundaryBehavior::Error, BoundaryBehavior::Error)
  }

  pub fn interp_with_extrap(
    &self,
    xi: T,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();

    if xi < x_min {
      if ulps_eq!(xi, x_min, max_ulps = 4) {
        return Ok(self.y[0]);
      }
      return self.extrapolate(left_extrap, xi, Side::Left);
    }

    if xi > x_max {
      let n = self.grid.n_points();
      if ulps_eq!(xi, x_max, max_ulps = 4) {
        return Ok(self.y[n - 1]);
      }
      return self.extrapolate(right_extrap, xi, Side::Right);
    }

    let idx = self.grid.find_interval_index(xi);
    Ok(self.interpolate_at(xi, idx))
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  fn edge(&self, side: Side) -> GridEdge
  where
    T: Float,
  {
    let (t, y) = match side {
      Side::Left => (self.grid.x_min(), self.y[0]),
      Side::Right => (self.grid.x_max(), self.y[self.grid.n_points() - 1]),
    };
    GridEdge {
      t: t.to_f64().unwrap(),
      y: y.to_f64().unwrap(),
    }
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  fn extrapolate(&self, behavior: BoundaryBehavior, xi: T, side: Side) -> Result<T, Report>
  where
    T: Float,
  {
    match behavior {
      BoundaryBehavior::Linear(law) => {
        let value = law.eval(self.edge(side), xi.to_f64().unwrap());
        Ok(T::from(value).unwrap())
      },
      BoundaryBehavior::HardApproach(law) => {
        let xi_f64 = xi.to_f64().unwrap();
        if (side == Side::Left && xi_f64 < law.t_hard) || (side == Side::Right && xi_f64 > law.t_hard) {
          return Ok(T::zero());
        }
        let value = law.eval(self.edge(side), xi_f64);
        Ok(T::from(value).unwrap())
      },
      BoundaryBehavior::Hard => Ok(T::zero()),
      BoundaryBehavior::Error => {
        let (side_word, bound) = match side {
          Side::Left => ("below", self.grid.x_min()),
          Side::Right => ("above", self.grid.x_max()),
        };
        make_error!(
          "GridFn evaluated at {xi:?}, {side_word} the support boundary {bound:?}, but no extrapolation policy is set for that side"
        )
      },
    }
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
  pub fn shift_y(&self, delta: T) -> Self
  where
    T: Float,
  {
    Self {
      grid: self.grid,
      y: self.y.mapv(|v| v + delta),
    }
  }

  #[allow(clippy::integer_division, reason = "integer division is the intended floor division")]
  pub fn negate_arg_inplace(&mut self) -> Result<(), Report>
  where
    T: Float,
  {
    let x_max = self.grid.x_max();
    let dx = self.grid.dx();
    let n_points = self.grid.n_points();
    self.grid = Grid::from_start_dx(-x_max, dx, n_points)?;

    let n = self.y.len();
    for i in 0..n / 2 {
      self.y.swap(i, n - 1 - i);
    }

    Ok(())
  }

  fn resample(&self, grid: &Grid<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    self.resample_with_extrap(grid, BoundaryBehavior::Error, BoundaryBehavior::Error)
  }

  pub fn resample_with_extrap(
    &self,
    grid: &Grid<T>,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let n_points = grid.n_points();
    let y_new = (0..n_points)
      .map(|i| self.interp_with_extrap(grid.x_at(i), left_extrap, right_extrap))
      .collect::<Result<Vec<T>, Report>>()?;
    Self::from_grid_array(*grid, Array1::from_vec(y_new))
  }

  pub fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let (x_min, x_max) = x_range;
    let grid = Grid::from_range_n_points(x_min, x_max, n_points)?;
    self.resample(&grid)
  }

  pub fn resample_range_dx_clamped(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_dx(x_range.0, x_range.1, dx)?;
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();
    let y_new = (0..grid.n_points())
      .map(|i| self.interp(grid.x_at(i).max(x_min).min(x_max)))
      .collect::<Result<Vec<T>, Report>>()?;
    Self::from_grid_array(grid, Array1::from_vec(y_new))
  }
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
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
