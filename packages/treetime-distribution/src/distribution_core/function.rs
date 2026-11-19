use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use ndarray_stats::errors::MinMaxError;
use treetime_grid::grid::Grid;
use treetime_grid::{BoundaryBehavior, EdgeCell, GridFn, HardCell, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::float_compare::within_epsilons_of_scale;
use treetime_utils::make_error;

const MAX_COST_ULPS: f64 = 10.0;

#[derive(Clone, Debug)]
pub struct DistributionFunction {
  grid_fn: GridFn,
  ampl: f64,
  left_extrap: BoundaryBehavior,
  right_extrap: BoundaryBehavior,
}

impl DistributionFunction {
  pub(crate) fn from_arrays(x: &Array1<f64>, y: Array1<f64>) -> Result<Self, Report> {
    let grid_fn = GridFn::from_arrays(x, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_range_values(x_range: (f64, f64), y: Array1<f64>) -> Result<Self, Report> {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  #[cfg_attr(
    dylint_lib = "custom",
    allow(
      pub_unused_in_workspace,
      reason = "used only by tests of other workspace crates, which a cfg(test) item cannot reach"
    )
  )]
  pub fn from_start_dx_values(x_min: f64, dx: f64, y: Array1<f64>) -> Result<Self, Report> {
    let grid_fn = GridFn::from_start_dx_values(x_min, dx, y)?;
    Ok(Self::from_grid_fn(grid_fn))
  }

  pub fn from_grid_fn(grid_fn: GridFn) -> Self {
    Self {
      grid_fn,
      ampl: 0.0,
      left_extrap: BoundaryBehavior::default(),
      right_extrap: BoundaryBehavior::default(),
    }
  }

  fn from_grid_fn_with_extrap(
    grid_fn: GridFn,
    ampl: f64,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Self {
    Self {
      grid_fn,
      ampl,
      left_extrap,
      right_extrap,
    }
  }

  pub(crate) fn grid_fn(&self) -> &GridFn {
    &self.grid_fn
  }

  pub fn t(&self) -> Array1<f64> {
    self.grid_fn.x()
  }

  pub fn x_min(&self) -> f64 {
    self.grid_fn.x_min()
  }

  pub fn x_max(&self) -> f64 {
    self.grid_fn.x_max()
  }

  pub fn dx(&self) -> f64 {
    self.grid_fn.dx()
  }

  pub fn x_at(&self, idx: usize) -> f64 {
    self.grid_fn.x_at(idx)
  }

  pub fn y(&self) -> &Array1<f64> {
    self.grid_fn.y()
  }

  pub fn amplitude(&self) -> f64 {
    self.ampl
  }

  #[must_use]
  pub fn with_amplitude(mut self, ampl: f64) -> Self {
    self.ampl = ampl;
    self
  }

  pub fn values(&self) -> Array1<f64> {
    self.y().mapv(|v| self.ampl + v)
  }

  pub(crate) fn interp(&self, x: f64) -> Result<f64, Report> {
    Ok(self.ampl + self.interp_stored(x)?)
  }

  pub(crate) fn interp_many(&self, xs: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let values = xs
      .iter()
      .map(|&q| self.interp(q))
      .collect::<Result<Vec<f64>, Report>>()?;
    Ok(Array1::from_vec(values))
  }

  pub(crate) fn interp_stored(&self, x: f64) -> Result<f64, Report> {
    let val = self
      .grid_fn
      .interp_with_extrap(x, self.left_extrap, self.right_extrap)?;
    if val == 0.0 && self.is_beyond_hard_boundary(x) {
      return Ok(f64::INFINITY);
    }
    Ok(val)
  }

  pub(crate) fn interp_many_stored(&self, xs: &Array1<f64>) -> Result<Array1<f64>, Report> {
    let values = xs
      .iter()
      .map(|&q| self.interp_stored(q))
      .collect::<Result<Vec<f64>, Report>>()?;
    Ok(Array1::from_vec(values))
  }

  fn is_beyond_hard_boundary(&self, x: f64) -> bool {
    let x_min = self.x_min();
    let x_max = self.x_max();
    if x < x_min {
      return match self.left_extrap() {
        BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_) => true,
        BoundaryBehavior::HardApproach(law) => x < law.t_hard,
        _ => false,
      };
    }
    if x > x_max {
      return match self.right_extrap() {
        BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_) => true,
        BoundaryBehavior::HardApproach(law) => x > law.t_hard,
        _ => false,
      };
    }
    false
  }

  pub fn left_extrap(&self) -> BoundaryBehavior {
    self.left_extrap
  }

  pub fn right_extrap(&self) -> BoundaryBehavior {
    self.right_extrap
  }

  pub fn extrap(&self, side: Side) -> BoundaryBehavior {
    match side {
      Side::Left => self.left_extrap,
      Side::Right => self.right_extrap,
    }
  }

  pub fn hard_cell(&self, cell: usize) -> Result<Option<HardCell>, Report> {
    self.grid_fn.hard_cell(cell, self.left_extrap, self.right_extrap)
  }

  pub fn edge_cell(&self, cell: usize) -> Result<Option<EdgeCell>, Report> {
    self.grid_fn.edge_cell(cell, self.left_extrap, self.right_extrap)
  }

  pub fn with_left_extrap(mut self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    self.left_extrap = behavior;
    Ok(self)
  }

  pub fn with_right_extrap(mut self, behavior: BoundaryBehavior) -> Result<Self, Report> {
    self.right_extrap = behavior;
    Ok(self)
  }

  fn resample(&self, grid: &Grid) -> Result<Self, Report> {
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

  pub(crate) fn resample_range_n_points(&self, x_range: (f64, f64), n_points: usize) -> Result<Self, Report> {
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
  pub fn resample_range_dx(&self, x_range: (f64, f64), dx: f64, max_points: MaxGridPoints) -> Result<Self, Report> {
    let grid = Grid::from_range_dx(x_range.0, x_range.1, dx, max_points)?;
    self.resample(&grid)
  }

  fn resample_range_dx_clamped(
    &self,
    x_range: (f64, f64),
    dx: f64,
    anchor: Side,
    max_points: MaxGridPoints,
  ) -> Result<Self, Report> {
    let grid_fn = self
      .grid_fn
      .resample_range_dx_clamped(x_range, dx, anchor, max_points)?;
    Ok(Self::from_grid_fn_with_extrap(
      grid_fn,
      self.ampl,
      self.left_extrap,
      self.right_extrap,
    ))
  }

  pub(crate) fn resample_dx(&self, dx: f64, max_points: MaxGridPoints) -> Result<Self, Report> {
    let y = self.y();
    let anchor = resample_anchor(self.left_extrap, self.right_extrap, y[0], y[y.len() - 1]);
    self.resample_range_dx_clamped((self.x_min(), self.x_max()), dx, anchor, max_points)
  }

  pub(crate) fn len(&self) -> usize {
    self.grid_fn.len()
  }

  pub(crate) fn is_empty(&self) -> bool {
    self.grid_fn.is_empty()
  }

  pub(crate) fn negate_arg_inplace(&mut self) -> Result<(), Report> {
    self.grid_fn.negate_arg_inplace()?;
    self.left_extrap = negate_tail_law(self.left_extrap);
    self.right_extrap = negate_tail_law(self.right_extrap);
    std::mem::swap(&mut self.left_extrap, &mut self.right_extrap);
    Ok(())
  }

  pub(crate) fn likely_time(&self) -> Result<Option<f64>, Report> {
    match self.y().argmin() {
      Ok(idx) => Ok(Some(self.t()[idx])),
      Err(MinMaxError::EmptyInput) => Ok(None),
      Err(MinMaxError::UndefinedOrder) => {
        make_error!("Cannot find the most likely time of a distribution function: its values contain NaN")
      },
    }
  }

  #[must_use]
  pub(crate) fn shift_y(&self, delta: f64) -> Self {
    Self::from_grid_fn_with_extrap(
      self.grid_fn.shift_y(delta),
      self.ampl,
      self.left_extrap,
      self.right_extrap,
    )
  }

  pub fn scale_time(&self, factor: f64) -> Result<Self, Report> {
    Ok(Self::from_grid_fn_with_extrap(
      self.grid_fn.scale_x(factor)?,
      self.ampl,
      scale_tail_law(self.left_extrap, factor),
      scale_tail_law(self.right_extrap, factor),
    ))
  }

  pub fn shift_time(&self, shift: f64) -> Result<Self, Report> {
    Ok(Self::from_grid_fn_with_extrap(
      self.grid_fn.shift_x(shift)?,
      self.ampl,
      shift_tail_law(self.left_extrap, shift),
      shift_tail_law(self.right_extrap, shift),
    ))
  }

  pub fn repeak(&self) -> Result<Option<Self>, Report> {
    let peak = match self.y().min() {
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
    if !peak.is_finite() {
      return Ok(None);
    }
    let grid_fn = self.grid_fn.map_y(|v| v - peak);
    Ok(Some(Self::from_grid_fn_with_extrap(
      grid_fn,
      self.ampl + peak,
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

impl PartialEq for DistributionFunction {
  fn eq(&self, other: &Self) -> bool {
    let with_level = |f: &Self| f.grid_fn.map_y(|v| f.ampl + v);
    self.left_extrap == other.left_extrap
      && self.right_extrap == other.right_extrap
      && with_level(self) == with_level(other)
  }
}

fn resample_anchor(left: BoundaryBehavior, right: BoundaryBehavior, first: f64, last: f64) -> Side {
  let right_end_denser = last < first;
  match (left.is_soft(), right.is_soft()) {
    (false, true) => Side::Left,
    (true, false) => Side::Right,
    _ if right_end_denser && !costs_coincide(first, last) => Side::Right,
    _ => Side::Left,
  }
}

pub(crate) fn costs_coincide(a: f64, b: f64) -> bool {
  within_epsilons_of_scale(a, b, a.abs().max(b.abs()).max(1.0), MAX_COST_ULPS)
}

fn shift_tail_law(behavior: BoundaryBehavior, shift: f64) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.shift_arg(shift)),
    BoundaryBehavior::HardEdgeCell(law) => BoundaryBehavior::HardEdgeCell(law.shift_arg(shift)),
    BoundaryBehavior::Error | BoundaryBehavior::Hard | BoundaryBehavior::Linear(_) => behavior,
  }
}

fn negate_tail_law(behavior: BoundaryBehavior) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.negate_arg()),
    BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(law.negate_arg()),
    BoundaryBehavior::HardEdgeCell(law) => BoundaryBehavior::HardEdgeCell(law.negate_arg()),
    BoundaryBehavior::Error | BoundaryBehavior::Hard => behavior,
  }
}

fn scale_tail_law(behavior: BoundaryBehavior, factor: f64) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(law.scale_arg(factor)),
    BoundaryBehavior::Linear(law) => BoundaryBehavior::Linear(law.scale_arg(factor)),
    BoundaryBehavior::HardEdgeCell(law) => BoundaryBehavior::HardEdgeCell(law.scale_arg(factor)),
    other => other,
  }
}
