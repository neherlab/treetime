use crate::policy::{Plain, PolicyMarker, YAxisPolicy};
use approx::UlpsEq;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use ndarray_stats::errors::MinMaxError;
use num::Float;
use treetime_grid::grid::Grid;
use treetime_grid::{BoundaryBehavior, GridFn, InterpElem, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::make_error;

#[derive(Clone, Debug)]
pub struct DistributionFunction<T: InterpElem, Y: YAxisPolicy = Plain> {
  grid_fn: GridFn<T>,
  ampl: T,
  left_extrap: BoundaryBehavior,
  right_extrap: BoundaryBehavior,
  _policy: PolicyMarker<Y>,
}

impl<T: InterpElem, Y: YAxisPolicy> PartialEq for DistributionFunction<T, Y> {
  fn eq(&self, other: &Self) -> bool {
    let with_level = |f: &Self| {
      let ampl: f64 = f.ampl.into();
      f.grid_fn.map_y(|v| Y::multiply(ampl, v.into()).into())
    };
    self.left_extrap == other.left_extrap
      && self.right_extrap == other.right_extrap
      && with_level(self) == with_level(other)
  }
}

impl<T: InterpElem, Y: YAxisPolicy> DistributionFunction<T, Y> {
  pub(crate) fn from_arrays(x: &Array1<T>, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays(x, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub(crate) fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>, max_points: MaxGridPoints) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays_nonuniform(x, y, max_points)?;
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

  pub fn from_grid_fn(grid_fn: GridFn<T>) -> Self {
    Self {
      grid_fn,
      ampl: unit_amplitude::<T, Y>(),
      left_extrap: BoundaryBehavior::default(),
      right_extrap: BoundaryBehavior::default(),
      _policy: PolicyMarker::new(),
    }
  }

  fn from_grid_fn_with_extrap(
    grid_fn: GridFn<T>,
    ampl: T,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Self {
    Self {
      grid_fn,
      ampl,
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

  pub fn amplitude(&self) -> T {
    self.ampl
  }

  #[must_use]
  pub fn with_amplitude(mut self, ampl: T) -> Self {
    self.ampl = ampl;
    self
  }

  pub fn values(&self) -> Array1<T>
  where
    T: Float,
  {
    let ampl: f64 = self.ampl.into();
    self.y().mapv(|v| Y::multiply(ampl, v.into()).into())
  }

  pub(crate) fn interp(&self, x: T) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    Ok(Y::multiply(self.ampl.into(), self.interp_stored(x)?.into()).into())
  }

  pub(crate) fn interp_many(&self, xs: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float + UlpsEq,
  {
    let values = xs.iter().map(|&q| self.interp(q)).collect::<Result<Vec<T>, Report>>()?;
    Ok(Array1::from_vec(values))
  }

  pub(crate) fn interp_stored(&self, x: T) -> Result<T, Report>
  where
    T: Float + UlpsEq,
  {
    let val = self
      .grid_fn
      .interp_with_extrap(x, self.left_extrap, self.right_extrap)?;
    let prob_zero: T = Y::probability_zero().into();
    if prob_zero != T::zero() && val == T::zero() && self.is_beyond_hard_boundary(x) {
      return Ok(prob_zero);
    }
    Ok(val)
  }

  pub(crate) fn interp_many_stored(&self, xs: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float + UlpsEq,
  {
    let values = xs
      .iter()
      .map(|&q| self.interp_stored(q))
      .collect::<Result<Vec<T>, Report>>()?;
    Ok(Array1::from_vec(values))
  }

  fn is_beyond_hard_boundary(&self, x: T) -> bool
  where
    T: Float,
  {
    let x_f64: f64 = x.into();
    let x_min: f64 = self.x_min().into();
    let x_max: f64 = self.x_max().into();
    if x_f64 < x_min {
      return match self.left_extrap() {
        BoundaryBehavior::Hard => true,
        BoundaryBehavior::HardApproach(law) => x_f64 < law.t_hard,
        _ => false,
      };
    }
    if x_f64 > x_max {
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
      self.ampl,
      self.left_extrap,
      self.right_extrap,
    ))
  }

  pub(crate) fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_n_points(x_range.0, x_range.1, n_points)?;
    self.resample(&grid)
  }

  #[cfg_attr(
    dylint_lib = "custom",
    allow(
      pub_unused_in_workspace,
      reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach"
    )
  )]
  pub fn resample_range_dx(&self, x_range: (T, T), dx: T, max_points: MaxGridPoints) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid = Grid::from_range_dx(x_range.0, x_range.1, dx, max_points)?;
    self.resample(&grid)
  }

  fn resample_range_dx_clamped(&self, x_range: (T, T), dx: T, max_points: MaxGridPoints) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_range_dx_clamped(x_range, dx, max_points)?;
    Ok(Self::from_grid_fn_with_extrap(
      grid_fn,
      self.ampl,
      self.left_extrap,
      self.right_extrap,
    ))
  }

  pub(crate) fn resample_dx(&self, dx: T, max_points: MaxGridPoints) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    self.resample_range_dx_clamped((self.x_min(), self.x_max()), dx, max_points)
  }

  pub(crate) fn len(&self) -> usize {
    self.grid_fn.len()
  }

  pub(crate) fn is_empty(&self) -> bool {
    self.grid_fn.is_empty()
  }

  pub(crate) fn negate_arg_inplace(&mut self) -> Result<(), Report>
  where
    T: Float,
  {
    self.grid_fn.negate_arg_inplace()?;
    self.left_extrap = scale_tail_law(self.left_extrap, -1.0);
    self.right_extrap = scale_tail_law(self.right_extrap, -1.0);
    std::mem::swap(&mut self.left_extrap, &mut self.right_extrap);
    Ok(())
  }

  pub(crate) fn likely_time(&self) -> Result<Option<T>, Report>
  where
    T: Float,
  {
    let y_values = self.y();
    let extremum = if Y::likely_is_maximum() {
      y_values.argmax()
    } else {
      y_values.argmin()
    };
    match extremum {
      Ok(idx) => Ok(Some(self.t()[idx])),
      Err(MinMaxError::EmptyInput) => Ok(None),
      Err(MinMaxError::UndefinedOrder) => {
        make_error!("Cannot find the most likely time of a distribution function: its values contain NaN")
      },
    }
  }

  #[must_use]
  pub(crate) fn shift_y(&self, delta: T) -> Self
  where
    T: Float,
  {
    Self::from_grid_fn_with_extrap(
      self.grid_fn.shift_y(delta),
      self.ampl,
      self.left_extrap,
      self.right_extrap,
    )
  }
}

impl<Y: YAxisPolicy> DistributionFunction<f64, Y> {
  pub fn scale_time(&self, factor: f64) -> Result<Self, Report> {
    Ok(Self::from_grid_fn_with_extrap(
      self.grid_fn.scale_x(factor)?,
      self.ampl,
      scale_tail_law(self.left_extrap, factor),
      scale_tail_law(self.right_extrap, factor),
    ))
  }

  pub fn repeak(&self) -> Result<Option<Self>, Report> {
    let stored = self.y();
    let extremum = if Y::likely_is_maximum() {
      stored.max()
    } else {
      stored.min()
    };
    let peak = match extremum {
      Ok(&peak) => peak,
      Err(MinMaxError::EmptyInput) => return Ok(None),
      Err(MinMaxError::UndefinedOrder) => {
        return make_error!(
          "Cannot re-peak a distribution on [{}, {}]: its values contain NaN",
          self.x_min(),
          self.x_max()
        );
      },
    };
    if !Y::is_defined(peak) || !peak.is_finite() {
      return Ok(None);
    }
    let grid_fn = self.grid_fn.map_y(|v| Y::divide(v, peak));
    Ok(Some(Self::from_grid_fn_with_extrap(
      grid_fn,
      Y::multiply(self.ampl, peak),
      self.left_extrap,
      self.right_extrap,
    )))
  }

  pub(crate) fn fit_soft_tail(self, side: Side, n_fit: usize) -> Result<Self, Report> {
    let law = SoftTailLaw::fit(&self.grid_fn, side, n_fit)?;
    match side {
      Side::Left => self.with_left_extrap(BoundaryBehavior::Linear(law)),
      Side::Right => self.with_right_extrap(BoundaryBehavior::Linear(law)),
    }
  }
}

fn unit_amplitude<T: InterpElem, Y: YAxisPolicy>() -> T {
  Y::from_plain(1.0).into()
}

fn scale_tail_law(behavior: BoundaryBehavior, factor: f64) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.scale_arg(factor)),
    BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(law.scale_arg(factor)),
    other => other,
  }
}
