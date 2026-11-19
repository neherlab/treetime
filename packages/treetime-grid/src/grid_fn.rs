use crate::InterpElem;
use crate::boundary_behavior::BoundaryBehavior;
use crate::boundary_behavior::DEFAULT_TAIL_FIT_POINTS;
use crate::edge_cell::EdgeCell;
use crate::grid::Grid;
use crate::grid_edge::GridEdge;
use crate::hard_approach_law::{HardApproachLaw, Side};
use crate::hard_cell::HardCell;
use crate::max_grid_points::MaxGridPoints;
use crate::y_axis_policy::YAxisPolicy;
use eyre::Report;
use ndarray::{Array1, s};
use num::Float;
use std::fmt::Debug;
use std::ops::RangeInclusive;
use treetime_utils::make_error;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GridFn<T: InterpElem> {
  grid: Grid<T>,
  y: Array1<T>,
}

impl<T: InterpElem> GridFn<T> {
  pub fn from_grid_values(grid: Grid<T>, y: Array1<T>) -> Result<Self, Report> {
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
    T: Float,
  {
    if x.len() != y.len() {
      return make_error!(
        "x and y arrays must have the same length, got {} and {}",
        x.len(),
        y.len()
      );
    }
    let grid = Grid::from_array(x)?;
    Self::from_grid_values(grid, y)
  }

  pub fn from_start_dx_values(x_min: T, dx: T, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_start_dx(x_min, dx, y.len())?;
    Self::from_grid_values(grid, y)
  }

  pub fn from_range_values((x_min, x_max): (T, T), y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_range_n_points(x_min, x_max, y.len())?;
    Self::from_grid_values(grid, y)
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

  pub fn with_values(&self, y: Array1<T>) -> Result<Self, Report> {
    Self::from_grid_values(self.grid, y)
  }

  pub fn slice(&self, points: RangeInclusive<usize>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = self.grid.slice(points.clone())?;
    Self::from_grid_values(grid, self.y.slice(s![points]).to_owned())
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

  pub fn x_at(&self, idx: usize) -> T
  where
    T: Float,
  {
    self.grid.x_at(idx)
  }

  pub fn coincides(&self, a: f64, b: f64) -> bool
  where
    T: Float,
  {
    self.grid.coincides(a, b)
  }

  pub fn y_min(&self) -> T {
    self
      .y
      .iter()
      .copied()
      .skip(1)
      .fold(self.y[0], |a, b| if a < b { a } else { b })
  }

  pub fn interp<Y: YAxisPolicy>(&self, xi: T) -> Result<T, Report>
  where
    T: Float,
  {
    self.interp_with_extrap::<Y>(xi, BoundaryBehavior::Error, BoundaryBehavior::Error)
  }

  pub fn interp_with_extrap<Y: YAxisPolicy>(
    &self,
    xi: T,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Result<T, Report>
  where
    T: Float,
  {
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();

    if xi < x_min {
      if self.grid.coincides(xi.into(), x_min.into()) {
        return Ok(self.y[0]);
      }
      return self.extrapolate::<Y>(left_extrap, xi, Side::Left);
    }

    if xi > x_max {
      let n = self.grid.n_points();
      if self.grid.coincides(xi.into(), x_max.into()) {
        return Ok(self.y[n - 1]);
      }
      return self.extrapolate::<Y>(right_extrap, xi, Side::Right);
    }

    let idx = self.grid.find_interval_index(xi);
    self.interpolate_at::<Y>(xi, idx, left_extrap, right_extrap)
  }

  pub fn hard_cell<Y: YAxisPolicy>(
    &self,
    cell: usize,
    left: BoundaryBehavior,
    right: BoundaryBehavior,
  ) -> Result<Option<HardCell>, Report>
  where
    T: Float,
  {
    let n = self.grid.n_points();
    let Some(HardCellEnds { boundary, anchor, side }) = self.hard_cell_ends::<Y>(cell)? else {
      return Ok(None);
    };
    let stored = match side {
      Side::Left => (cell == 0).then_some(left),
      Side::Right => (cell + 1 == n - 1).then_some(right),
    };
    let t_boundary: f64 = self.grid.x_at(boundary).into();
    let law = match stored {
      Some(BoundaryBehavior::HardApproach(law)) if self.grid.coincides(law.t_hard, t_boundary) => Ok(law),
      _ => self.fit_hard_cell_law::<Y>(t_boundary, side, anchor),
    }?;
    let anchor = self.neg_log_point::<Y>(anchor);
    Ok(Some(HardCell { law, anchor }))
  }

  pub fn edge_cell<Y: YAxisPolicy>(
    &self,
    cell: usize,
    left: BoundaryBehavior,
    right: BoundaryBehavior,
  ) -> Result<Option<EdgeCell>, Report>
  where
    T: Float,
  {
    let n = self.grid.n_points();
    if cell + 1 >= n {
      return make_error!("Cell {cell} is outside the grid of {n} points");
    }
    let (law, edge, anchor) = match (left, right) {
      (BoundaryBehavior::HardEdgeCell(law), _) if cell == 0 => (law, 0, 1),
      (_, BoundaryBehavior::HardEdgeCell(law)) if cell + 2 == n => (law, n - 1, n - 2),
      _ => return Ok(None),
    };
    let (edge, anchor) = (self.neg_log_point::<Y>(edge), self.neg_log_point::<Y>(anchor));
    if !(self.grid.coincides(law.t_edge, edge.t) && edge.y.is_finite() && anchor.y.is_finite()) {
      return Ok(None);
    }
    EdgeCell::new(law, edge, anchor).map(Some)
  }

  pub fn hard_cell_fit_range<Y: YAxisPolicy>(&self, cell: usize) -> Result<Option<RangeInclusive<usize>>, Report>
  where
    T: Float,
  {
    Ok(
      self
        .hard_cell_ends::<Y>(cell)?
        .map(|HardCellEnds { anchor, side, .. }| self.fit_range::<Y>(anchor, side)),
    )
  }

  fn hard_cell_ends<Y: YAxisPolicy>(&self, cell: usize) -> Result<Option<HardCellEnds>, Report>
  where
    T: Float,
  {
    let n = self.grid.n_points();
    if cell + 1 >= n {
      return make_error!("Cell {cell} is outside the grid of {n} points");
    }
    Ok(match (self.is_hard_end::<Y>(cell), self.is_hard_end::<Y>(cell + 1)) {
      (true, false) => Some(HardCellEnds {
        boundary: cell,
        anchor: cell + 1,
        side: Side::Left,
      }),
      (false, true) => Some(HardCellEnds {
        boundary: cell + 1,
        anchor: cell,
        side: Side::Right,
      }),
      (true, true) | (false, false) => None,
    })
  }

  fn fit_range<Y: YAxisPolicy>(&self, anchor: usize, side: Side) -> RangeInclusive<usize>
  where
    T: Float,
  {
    let reach = DEFAULT_TAIL_FIT_POINTS - 1;
    let finite = |i: &usize| self.neg_log_at::<Y>(*i).is_finite();
    match side {
      Side::Right => {
        let start = (anchor.saturating_sub(reach)..=anchor)
          .rev()
          .take_while(finite)
          .last()
          .unwrap_or(anchor);
        start..=anchor
      },
      Side::Left => {
        let end = (anchor..=(anchor + reach).min(self.grid.n_points() - 1))
          .take_while(finite)
          .last()
          .unwrap_or(anchor);
        anchor..=end
      },
    }
  }

  fn neg_log_edge<Y: YAxisPolicy>(&self, side: Side) -> GridEdge
  where
    T: Float,
  {
    self.neg_log_point::<Y>(match side {
      Side::Left => 0,
      Side::Right => self.grid.n_points() - 1,
    })
  }

  fn neg_log_point<Y: YAxisPolicy>(&self, idx: usize) -> GridEdge
  where
    T: Float,
  {
    GridEdge {
      t: self.grid.x_at(idx).into(),
      y: self.neg_log_at::<Y>(idx),
    }
  }

  fn is_hard_end<Y: YAxisPolicy>(&self, idx: usize) -> bool {
    self.neg_log_at::<Y>(idx).is_infinite()
  }

  fn neg_log_at<Y: YAxisPolicy>(&self, idx: usize) -> f64 {
    Y::to_neg_log(self.y[idx].into())
  }

  fn extrapolate<Y: YAxisPolicy>(&self, behavior: BoundaryBehavior, xi: T, side: Side) -> Result<T, Report>
  where
    T: Float,
  {
    match behavior {
      BoundaryBehavior::Linear(law) => {
        let value = Y::from_neg_log(law.eval(self.neg_log_edge::<Y>(side), xi.into()));
        Ok(value.into())
      },
      BoundaryBehavior::HardApproach(law) => {
        let xi_f64: f64 = xi.into();
        if (side == Side::Left && xi_f64 < law.t_hard) || (side == Side::Right && xi_f64 > law.t_hard) {
          return Ok(T::zero());
        }
        let value = Y::from_neg_log(law.eval(self.neg_log_edge::<Y>(side), xi_f64));
        Ok(value.into())
      },
      BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_) => Ok(T::zero()),
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

  fn interpolate_at<Y: YAxisPolicy>(
    &self,
    q: T,
    idx: usize,
    left: BoundaryBehavior,
    right: BoundaryBehavior,
  ) -> Result<T, Report>
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
    let t = (q - x0) / (self.grid.x_at(idx + 1) - x0);
    let t_f64: f64 = t.into();
    if t_f64 <= 0.0 {
      return Ok(y0);
    }
    if t_f64 >= 1.0 {
      return Ok(y1);
    }
    if let Some(cell) = self.edge_cell::<Y>(idx, left, right)? {
      return Ok(<T as From<f64>>::from(Y::from_neg_log(cell.eval(q.into()))));
    }
    if t_f64 <= NODE_ROUNDING {
      return Ok(y0);
    }
    if t_f64 >= 1.0 - NODE_ROUNDING {
      return Ok(y1);
    }
    match self.hard_cell::<Y>(idx, left, right)? {
      Some(cell) => Ok(<T as From<f64>>::from(Y::from_neg_log(cell.eval(q.into())))),
      None if self.is_hard_end::<Y>(idx) || self.is_hard_end::<Y>(idx + 1) => Ok(y0.max(y1)),
      None => Ok(y0 + t * (y1 - y0)),
    }
  }

  fn fit_hard_cell_law<Y: YAxisPolicy>(
    &self,
    t_boundary: f64,
    side: Side,
    anchor: usize,
  ) -> Result<HardApproachLaw, Report>
  where
    T: Float,
  {
    let range = self.fit_range::<Y>(anchor, side);
    let indices: Vec<usize> = match side {
      Side::Right => range.rev().collect(),
      Side::Left => range.collect(),
    };
    let points: Vec<(f64, f64)> = indices
      .into_iter()
      .map(|i| {
        let point = self.neg_log_point::<Y>(i);
        (point.t, point.y)
      })
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
      grid: self.grid.scaled(factor)?,
      y: self.y.clone(),
    })
  }

  pub fn shift_x(&self, delta: T) -> Result<Self, Report>
  where
    T: Float,
  {
    Ok(Self {
      grid: self.grid.shifted(delta)?,
      y: self.y.clone(),
    })
  }

  #[expect(clippy::integer_division, reason = "integer division is the intended floor division")]
  pub fn negate_arg_inplace(&mut self) -> Result<(), Report>
  where
    T: Float,
  {
    self.grid = self.grid.negated()?;
    let n = self.y.len();
    for i in 0..n / 2 {
      self.y.swap(i, n - 1 - i);
    }
    Ok(())
  }

  pub fn resample_with_extrap<Y: YAxisPolicy>(
    &self,
    grid: &Grid<T>,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Result<Self, Report>
  where
    T: Float,
  {
    let n_points = grid.n_points();
    let y_new = (0..n_points)
      .map(|i| self.interp_with_extrap::<Y>(grid.x_at(i), left_extrap, right_extrap))
      .collect::<Result<Vec<T>, Report>>()?;
    Self::from_grid_values(*grid, Array1::from_vec(y_new))
  }

  pub fn resample_range_dx_clamped<Y: YAxisPolicy>(
    &self,
    x_range: (T, T),
    dx: T,
    anchor: Side,
    max_points: MaxGridPoints,
  ) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid = Grid::from_range_dx_anchored(x_range.0, x_range.1, dx, anchor, max_points)?;
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();
    let y_new = (0..grid.n_points())
      .map(|i| self.interp::<Y>(grid.x_at(i).max(x_min).min(x_max)))
      .collect::<Result<Vec<T>, Report>>()?;
    Self::from_grid_values(grid, Array1::from_vec(y_new))
  }
}

struct HardCellEnds {
  boundary: usize,
  anchor: usize,
  side: Side,
}
