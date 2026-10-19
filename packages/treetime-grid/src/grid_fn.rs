use crate::InterpElem;
use crate::boundary_behavior::BoundaryBehavior;
use crate::boundary_behavior::DEFAULT_TAIL_FIT_POINTS;
use crate::grid::{Grid, from_count};
use crate::grid_edge::GridEdge;
use crate::hard_approach_law::{HardApproachLaw, Side};
use crate::interp_nonuniform::interp_nonuniform;
use crate::max_grid_points::MaxGridPoints;
use approx::{UlpsEq, ulps_eq};
use eyre::{Report, WrapErr};
use ndarray::{Array1, s};
use ndarray_stats::QuantileExt;
use num::Float;
use std::fmt::Debug;
use treetime_utils::array::ndarray::has_uniform_spacing;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GridFn<T: InterpElem> {
  grid: Grid<T>,
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

  pub fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>, max_points: MaxGridPoints) -> Result<Self, Report>
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

    let x_min = x[0];
    let x_max = x[x.len() - 1];
    let range = (x_min.into(), x_max.into());
    let input_dx = (range.1 - range.0) / (x.len() - 1).to_f64_exact();
    max_points.point_count(x.len().to_f64_exact(), range, input_dx)?;

    if has_uniform_spacing(x) {
      let grid = Grid::from_array(x)?;
      return Self::from_grid_array(grid, y.clone());
    }

    let dx = find_min_spacing(x)?;
    let required = ((x_max - x_min) / dx).ceil() + T::one();
    let n_points = max_points.point_count(required.into(), range, dx.into())?;
    let dx = (x_max - x_min) / from_count::<T>(n_points - 1);

    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    let y_uniform = interp_nonuniform(x, y, n_points, |i| grid.x_at(i))?;
    Self::from_grid_array(grid, y_uniform)
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

  pub(crate) fn grid(&self) -> &Grid<T> {
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
    self.interpolate_at(xi, idx, left_extrap, right_extrap)
  }

  fn edge(&self, side: Side) -> GridEdge
  where
    T: Float,
  {
    let (t, y) = match side {
      Side::Left => (self.grid.x_min(), self.y[0]),
      Side::Right => (self.grid.x_max(), self.y[self.grid.n_points() - 1]),
    };
    GridEdge {
      t: t.into(),
      y: y.into(),
    }
  }

  fn extrapolate(&self, behavior: BoundaryBehavior, xi: T, side: Side) -> Result<T, Report>
  where
    T: Float,
  {
    match behavior {
      BoundaryBehavior::Linear(law) => {
        let value = law.eval(self.edge(side), xi.into());
        Ok(value.into())
      },
      BoundaryBehavior::HardApproach(law) => {
        let xi_f64: f64 = xi.into();
        if (side == Side::Left && xi_f64 < law.t_hard) || (side == Side::Right && xi_f64 > law.t_hard) {
          return Ok(T::zero());
        }
        let value = law.eval(self.edge(side), xi_f64);
        Ok(value.into())
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

  fn interpolate_at(&self, q: T, idx: usize, left: BoundaryBehavior, right: BoundaryBehavior) -> Result<T, Report>
  where
    T: Float,
  {
    const NODE_ROUNDING: f64 = 1e-9;
    let n = self.grid.n_points();
    if idx >= n - 1 {
      return Ok(self.y[n - 1]);
    }
    let x0 = self.grid.x_at(idx);
    let y0 = self.y[idx];
    let y1 = self.y[idx + 1];
    let t = (q - x0) / self.grid.dx();
    let t_f64: f64 = t.into();
    if t_f64 <= NODE_ROUNDING {
      return Ok(y0);
    }
    if t_f64 >= 1.0 - NODE_ROUNDING {
      return Ok(y1);
    }
    match (y0.is_infinite(), y1.is_infinite()) {
      (false, false) => Ok(y0 + t * (y1 - y0)),
      (true, true) => Ok(y0.max(y1)),
      (true, false) => self.hard_cell_value(q, idx, idx + 1, left, right),
      (false, true) => self.hard_cell_value(q, idx + 1, idx, left, right),
    }
  }

  fn hard_cell_value(
    &self,
    q: T,
    boundary: usize,
    anchor: usize,
    left: BoundaryBehavior,
    right: BoundaryBehavior,
  ) -> Result<T, Report>
  where
    T: Float,
  {
    let n = self.grid.n_points();
    let t_boundary: f64 = self.grid.x_at(boundary).into();
    let (side, stored) = if boundary > anchor {
      (Side::Right, (boundary == n - 1).then_some(right))
    } else {
      (Side::Left, (boundary == 0).then_some(left))
    };
    let stored_law = match stored {
      Some(BoundaryBehavior::HardApproach(law)) if ulps_eq!(law.t_hard, t_boundary, max_ulps = 4) => Some(law),
      _ => None,
    };
    let law = match stored_law {
      Some(law) => law,
      None => self.fit_hard_cell_law(t_boundary, side, anchor)?,
    };
    let edge = GridEdge {
      t: self.grid.x_at(anchor).into(),
      y: self.y[anchor].into(),
    };
    Ok(<T as From<f64>>::from(law.eval(edge, q.into())))
  }

  fn fit_hard_cell_law(&self, t_boundary: f64, side: Side, anchor: usize) -> Result<HardApproachLaw, Report>
  where
    T: Float,
  {
    let n = self.grid.n_points();
    let points: Vec<(f64, f64)> = (0..DEFAULT_TAIL_FIT_POINTS)
      .map_while(|k| match side {
        Side::Right => anchor.checked_sub(k),
        Side::Left => Some(anchor + k).filter(|&i| i < n),
      })
      .map(|i| (self.grid.x_at(i).into(), self.y[i].into()))
      .take_while(|(_, y): &(f64, f64)| y.is_finite())
      .collect();
    if points.len() < 2 {
      return Ok(HardApproachLaw::flat(t_boundary));
    }
    HardApproachLaw::from_points(t_boundary, side, &points)
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

  #[must_use]
  pub fn map_y(&self, f: impl Fn(T) -> T) -> Self {
    Self {
      grid: self.grid,
      y: self.y.mapv(f),
    }
  }

  pub fn scale_x(&self, factor: T) -> Result<Self, Report>
  where
    T: Float,
  {
    if factor.is_nan() || factor <= T::zero() {
      return make_error!("Argument scale factor must be positive, got {factor:?}");
    }
    Ok(Self {
      grid: Grid::from_start_dx(
        self.grid.x_min() * factor,
        self.grid.dx() * factor,
        self.grid.n_points(),
      )?,
      y: self.y.clone(),
    })
  }

  #[expect(clippy::integer_division, reason = "integer division is the intended floor division")]
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

  pub fn resample_range_dx_clamped(&self, x_range: (T, T), dx: T, max_points: MaxGridPoints) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_dx(x_range.0, x_range.1, dx, max_points)?;
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();
    let y_new = (0..grid.n_points())
      .map(|i| self.interp(grid.x_at(i).max(x_min).min(x_max)))
      .collect::<Result<Vec<T>, Report>>()?;
    Self::from_grid_array(grid, Array1::from_vec(y_new))
  }
}

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

  let min_dx = *diffs
    .min()
    .wrap_err("When finding the smallest spacing of the input array")?;

  if !min_dx.is_finite() || min_dx <= T::zero() {
    return make_error!("Invalid spacing in input array: {min_dx:?}");
  }

  Ok(min_dx)
}
