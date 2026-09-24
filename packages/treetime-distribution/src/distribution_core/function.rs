use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use approx::UlpsEq;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_grid::grid::Grid;
use treetime_grid::{BoundaryBehavior, GridFn, InterpElem, Side, SoftTailLaw};

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct DistributionFunction<T: InterpElem, Y: YAxisPolicy = Plain> {
  grid_fn: GridFn<T>,
  #[serde(skip)]
  left_extrap: BoundaryBehavior,
  #[serde(skip)]
  right_extrap: BoundaryBehavior,
  #[serde(skip)]
  _policy: PolicyMarker<Y>,
}

impl<T: InterpElem, Y: YAxisPolicy> DistributionFunction<T, Y> {
  pub(crate) fn from_arrays(x: &Array1<T>, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays(x, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub(crate) fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays_nonuniform(x, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_range_values(x_range: (T, T), y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_start_dx_values(x_min: T, dx: T, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::from_start_dx_values(x_min, dx, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_n_points<F>((x_min, x_max): (T, T), n_points: usize, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid_fn = GridFn::from_n_points((x_min, x_max), n_points, y_fn)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_grid<F>((x_min, x_max): (T, T), dx: T, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid_fn = GridFn::from_grid((x_min, x_max), dx, y_fn)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn constant((x_min, x_max): (T, T), n_points: usize, value: T) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::constant((x_min, x_max), n_points, value)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn zeros((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::zeros((x_min, x_max), n_points)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn ones((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::ones((x_min, x_max), n_points)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_grid_fn(grid_fn: GridFn<T>) -> Self {
    Self {
      grid_fn,
      left_extrap: BoundaryBehavior::default(),
      right_extrap: BoundaryBehavior::default(),
      _policy: PolicyMarker::new(),
    }
  }

  fn from_grid_fn_with_extrap(
    grid_fn: GridFn<T>,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Self {
    Self {
      grid_fn,
      left_extrap,
      right_extrap,
      _policy: PolicyMarker::new(),
    }
  }

  pub(crate) fn grid_fn(&self) -> &GridFn<T> {
    &self.grid_fn
  }

  pub(crate) fn t(&self) -> Array1<T>
  where
    T: Float,
  {
    self.grid_fn.x()
  }

  pub(crate) fn x_min(&self) -> T {
    self.grid_fn.x_min()
  }

  pub(crate) fn x_max(&self) -> T
  where
    T: Float,
  {
    self.grid_fn.x_max()
  }

  pub fn dx(&self) -> T {
    self.grid_fn.dx()
  }

  pub fn y(&self) -> &Array1<T> {
    self.grid_fn.y()
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  pub(crate) fn interp(&self, x: T) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    let val = self
      .grid_fn
      .interp_with_extrap(x, self.left_extrap, self.right_extrap)?;
    let prob_zero = T::from(Y::probability_zero()).unwrap();
    if prob_zero != T::zero() && val == T::zero() && self.is_beyond_hard_boundary(x) {
      return Ok(prob_zero);
    }
    Ok(val)
  }

  pub(crate) fn interp_many(&self, xs: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float + UlpsEq,
  {
    let values = xs.iter().map(|&q| self.interp(q)).collect::<Result<Vec<T>, Report>>()?;
    Ok(Array1::from_vec(values))
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  fn is_beyond_hard_boundary(&self, x: T) -> bool
  where
    T: Float,
  {
    let x_f64 = x.to_f64().unwrap();
    if x_f64 < self.x_min().to_f64().unwrap() {
      return match self.left_extrap() {
        BoundaryBehavior::Hard => true,
        BoundaryBehavior::HardApproach(law) => x_f64 < law.t_hard,
        _ => false,
      };
    }
    if x_f64 > self.x_max().to_f64().unwrap() {
      return match self.right_extrap() {
        BoundaryBehavior::Hard => true,
        BoundaryBehavior::HardApproach(law) => x_f64 > law.t_hard,
        _ => false,
      };
    }
    false
  }

  pub(crate) fn left_extrap(&self) -> BoundaryBehavior {
    self.left_extrap
  }

  pub fn right_extrap(&self) -> BoundaryBehavior {
    self.right_extrap
  }

  pub fn with_left_extrap(mut self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    self.left_extrap = behavior;
    Ok(self)
  }

  pub fn with_right_extrap(mut self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    self.right_extrap = behavior;
    Ok(self)
  }

  fn resample(&self, grid: &Grid<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self
      .grid_fn
      .resample_with_extrap(grid, self.left_extrap, self.right_extrap)?;
    Ok(Self::from_grid_fn_with_extrap(
      grid_fn,
      self.left_extrap,
      self.right_extrap,
    ))
  }

  pub fn resample_start_dx(&self, x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_start_dx(x_min, dx, n_points)?;
    self.resample(&grid)
  }

  pub(crate) fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_n_points(x_range.0, x_range.1, n_points)?;
    self.resample(&grid)
  }

  pub fn resample_range_dx(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_dx(x_range.0, x_range.1, dx)?;
    self.resample(&grid)
  }

  fn resample_range_dx_clamped(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_range_dx_clamped(x_range, dx)?;
    Ok(Self::from_grid_fn_with_extrap(
      grid_fn,
      self.left_extrap,
      self.right_extrap,
    ))
  }

  pub(crate) fn resample_dx(&self, dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    self.resample_range_dx_clamped((self.x_min(), self.x_max()), dx)
  }

  pub(crate) fn len(&self) -> usize {
    self.grid_fn.len()
  }

  pub(crate) fn is_empty(&self) -> bool {
    self.grid_fn.len() == 0
  }

  pub(crate) fn negate_arg_inplace(&mut self) -> Result<(), Report>
  where
    T: Float,
  {
    self.grid_fn.negate_arg_inplace()?;
    self.left_extrap = negate_tail_law(self.left_extrap);
    self.right_extrap = negate_tail_law(self.right_extrap);
    std::mem::swap(&mut self.left_extrap, &mut self.right_extrap);
    Ok(())
  }

  pub(crate) fn likely_time(&self) -> Option<T>
  where
    T: Float,
  {
    let y_values = self.y();
    let extremum = if Y::likely_is_maximum() {
      y_values.argmax()
    } else {
      y_values.argmin()
    };
    extremum.ok().map(|idx| self.t()[idx])
  }

  #[allow(
    clippy::unwrap_used,
    reason = "unwrap on a value an upstream invariant guarantees is present"
  )]
  #[must_use]
  pub(crate) fn shift_y(&self, delta: T) -> Self
  where
    T: Float,
  {
    Self::from_grid_fn_with_extrap(self.grid_fn.shift_y(delta), self.left_extrap, self.right_extrap)
  }
}

impl<Y: YAxisPolicy> DistributionFunction<f64, Y> {
  pub(crate) fn fit_soft_tail(self, side: Side, n_fit: usize) -> Result<Self, Report> {
    let law = SoftTailLaw::fit(&self.grid_fn, side, n_fit)?;
    match side {
      Side::Left => self.with_left_extrap(BoundaryBehavior::Linear(law)),
      Side::Right => self.with_right_extrap(BoundaryBehavior::Linear(law)),
    }
  }
}

fn negate_tail_law(behavior: BoundaryBehavior) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.negate_arg()),
    BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(law.negate_arg()),
    other => other,
  }
}
