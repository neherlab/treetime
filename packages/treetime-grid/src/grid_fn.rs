use crate::boundary_behavior::BoundaryBehavior;
use crate::boundary_behavior::DEFAULT_TAIL_FIT_POINTS;
use crate::edge_cell::EdgeCell;
use crate::grid::Grid;
use crate::grid_edge::GridEdge;
use crate::hard_approach_law::{HardApproachLaw, Side};
use crate::hard_cell::HardCell;
use crate::max_grid_points::MaxGridPoints;
use eyre::Report;
use ndarray::{Array1, s};
use std::fmt::Debug;
use std::ops::RangeInclusive;
use treetime_utils::make_error;

#[derive(Debug, Clone, PartialEq)]
pub struct GridFn {
  grid: Grid,
  y: Array1<f64>,
}

impl GridFn {
  pub fn from_grid_values(grid: Grid, y: Array1<f64>) -> Result<Self, Report> {
    if grid.n_points() != y.len() {
      return make_error!(
        "Grid has {} points but y array has {} elements",
        grid.n_points(),
        y.len()
      );
    }
    Ok(Self { grid, y })
  }

  pub fn from_arrays(x: &Array1<f64>, y: Array1<f64>) -> Result<Self, Report> {
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

  pub fn from_range_values((x_min, x_max): (f64, f64), y: Array1<f64>) -> Result<Self, Report> {
    let grid = Grid::from_range_n_points(x_min, x_max, y.len())?;
    Self::from_grid_values(grid, y)
  }

  // TODO: inefficient. Try to remove this method
  pub fn x(&self) -> Array1<f64> {
    self.grid.to_array()
  }

  pub fn y(&self) -> &Array1<f64> {
    &self.y
  }

  pub fn grid(&self) -> &Grid {
    &self.grid
  }

  pub fn with_values(&self, y: Array1<f64>) -> Result<Self, Report> {
    Self::from_grid_values(self.grid, y)
  }

  pub fn slice(&self, points: RangeInclusive<usize>) -> Result<Self, Report> {
    let grid = self.grid.slice(points.clone())?;
    Self::from_grid_values(grid, self.y.slice(s![points]).to_owned())
  }

  pub fn x_min(&self) -> f64 {
    self.grid.x_min()
  }

  pub fn x_max(&self) -> f64 {
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

  pub fn dx(&self) -> f64 {
    self.grid.dx()
  }

  pub fn x_at(&self, idx: usize) -> f64 {
    self.grid.x_at(idx)
  }

  pub fn coincides(&self, a: f64, b: f64) -> bool {
    self.grid.coincides(a, b)
  }

  pub fn y_min(&self) -> f64 {
    self
      .y
      .iter()
      .copied()
      .skip(1)
      .fold(self.y[0], |a, b| if a < b { a } else { b })
  }

  pub fn interp_with_extrap(
    &self,
    xi: f64,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Result<f64, Report> {
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();

    if xi < x_min {
      if self.grid.coincides(xi, x_min) {
        return Ok(self.y[0]);
      }
      return self.extrapolate(left_extrap, xi, Side::Left);
    }

    if xi > x_max {
      let n = self.grid.n_points();
      if self.grid.coincides(xi, x_max) {
        return Ok(self.y[n - 1]);
      }
      return self.extrapolate(right_extrap, xi, Side::Right);
    }

    let idx = self.grid.find_interval_index(xi);
    self.interpolate_at(xi, idx, left_extrap, right_extrap)
  }

  pub fn hard_cell(
    &self,
    cell: usize,
    left: BoundaryBehavior,
    right: BoundaryBehavior,
  ) -> Result<Option<HardCell>, Report> {
    let n = self.grid.n_points();
    let Some(HardCellEnds { boundary, anchor, side }) = self.hard_cell_ends(cell)? else {
      return Ok(None);
    };
    let stored = match side {
      Side::Left => (cell == 0).then_some(left),
      Side::Right => (cell + 1 == n - 1).then_some(right),
    };
    let t_boundary = self.grid.x_at(boundary);
    let law = match stored {
      Some(BoundaryBehavior::HardApproach(law)) if self.grid.coincides(law.t_hard, t_boundary) => Ok(law),
      _ => self.fit_hard_cell_law(t_boundary, side, anchor),
    }?;
    let anchor = self.neg_log_point(anchor);
    Ok(Some(HardCell { law, anchor }))
  }

  pub fn edge_cell(
    &self,
    cell: usize,
    left: BoundaryBehavior,
    right: BoundaryBehavior,
  ) -> Result<Option<EdgeCell>, Report> {
    let n = self.grid.n_points();
    if cell + 1 >= n {
      return make_error!("Cell {cell} is outside the grid of {n} points");
    }
    let (law, edge, anchor) = match (left, right) {
      (BoundaryBehavior::HardEdgeCell(law), _) if cell == 0 => (law, 0, 1),
      (_, BoundaryBehavior::HardEdgeCell(law)) if cell + 2 == n => (law, n - 1, n - 2),
      _ => return Ok(None),
    };
    let (edge, anchor) = (self.neg_log_point(edge), self.neg_log_point(anchor));
    if !(self.grid.coincides(law.t_edge, edge.t) && edge.y.is_finite() && anchor.y.is_finite()) {
      return Ok(None);
    }
    EdgeCell::new(law, edge, anchor).map(Some)
  }

  pub fn hard_cell_fit_range(&self, cell: usize) -> Result<Option<RangeInclusive<usize>>, Report> {
    Ok(
      self
        .hard_cell_ends(cell)?
        .map(|HardCellEnds { anchor, side, .. }| self.fit_range(anchor, side)),
    )
  }

  fn hard_cell_ends(&self, cell: usize) -> Result<Option<HardCellEnds>, Report> {
    let n = self.grid.n_points();
    if cell + 1 >= n {
      return make_error!("Cell {cell} is outside the grid of {n} points");
    }
    Ok(match (self.is_hard_end(cell), self.is_hard_end(cell + 1)) {
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

  fn fit_range(&self, anchor: usize, side: Side) -> RangeInclusive<usize> {
    let reach = DEFAULT_TAIL_FIT_POINTS - 1;
    let finite = |i: &usize| self.y[*i].is_finite();
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

  fn neg_log_edge(&self, side: Side) -> GridEdge {
    self.neg_log_point(match side {
      Side::Left => 0,
      Side::Right => self.grid.n_points() - 1,
    })
  }

  fn neg_log_point(&self, idx: usize) -> GridEdge {
    GridEdge {
      t: self.grid.x_at(idx),
      y: self.y[idx],
    }
  }

  fn is_hard_end(&self, idx: usize) -> bool {
    self.y[idx].is_infinite()
  }

  fn extrapolate(&self, behavior: BoundaryBehavior, xi: f64, side: Side) -> Result<f64, Report> {
    match behavior {
      BoundaryBehavior::Linear(law) => Ok(law.eval(self.neg_log_edge(side), xi)),
      BoundaryBehavior::HardApproach(law) => {
        if (side == Side::Left && xi < law.t_hard) || (side == Side::Right && xi > law.t_hard) {
          return Ok(f64::INFINITY);
        }
        Ok(law.eval(self.neg_log_edge(side), xi))
      },
      BoundaryBehavior::Hard | BoundaryBehavior::HardEdgeCell(_) => Ok(f64::INFINITY),
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

  fn interpolate_at(&self, q: f64, idx: usize, left: BoundaryBehavior, right: BoundaryBehavior) -> Result<f64, Report> {
    const NODE_ROUNDING: f64 = 1e-9;
    let n = self.grid.n_points();
    if idx >= n - 1 {
      return Ok(self.y[n - 1]);
    }
    let x0 = self.grid.x_at(idx);
    let y0 = self.y[idx];
    let y1 = self.y[idx + 1];
    let t = (q - x0) / (self.grid.x_at(idx + 1) - x0);
    if t <= 0.0 {
      return Ok(y0);
    }
    if t >= 1.0 {
      return Ok(y1);
    }
    if let Some(cell) = self.edge_cell(idx, left, right)? {
      return Ok(cell.eval(q));
    }
    if t <= NODE_ROUNDING {
      return Ok(y0);
    }
    if t >= 1.0 - NODE_ROUNDING {
      return Ok(y1);
    }
    match self.hard_cell(idx, left, right)? {
      Some(cell) => Ok(cell.eval(q)),
      None if self.is_hard_end(idx) || self.is_hard_end(idx + 1) => Ok(y0.max(y1)),
      None => Ok(y0 + t * (y1 - y0)),
    }
  }

  fn fit_hard_cell_law(&self, t_boundary: f64, side: Side, anchor: usize) -> Result<HardApproachLaw, Report> {
    let range = self.fit_range(anchor, side);
    let indices: Vec<usize> = match side {
      Side::Right => range.rev().collect(),
      Side::Left => range.collect(),
    };
    let points: Vec<(f64, f64)> = indices
      .into_iter()
      .map(|i| {
        let point = self.neg_log_point(i);
        (point.t, point.y)
      })
      .collect();
    if points.len() < 2 {
      return Ok(HardApproachLaw::flat(t_boundary));
    }
    HardApproachLaw::from_points(t_boundary, side, &points)
  }

  #[must_use]
  pub fn shift_y(&self, delta: f64) -> Self {
    Self {
      grid: self.grid,
      y: self.y.mapv(|v| v + delta),
    }
  }

  #[must_use]
  pub fn map_y(&self, f: impl Fn(f64) -> f64) -> Self {
    Self {
      grid: self.grid,
      y: self.y.mapv(f),
    }
  }

  pub fn scale_x(&self, factor: f64) -> Result<Self, Report> {
    if factor.is_nan() || factor <= 0.0 {
      return make_error!("Argument scale factor must be positive, got {factor:?}");
    }
    Ok(Self {
      grid: self.grid.scaled(factor)?,
      y: self.y.clone(),
    })
  }

  pub fn shift_x(&self, delta: f64) -> Result<Self, Report> {
    Ok(Self {
      grid: self.grid.shifted(delta)?,
      y: self.y.clone(),
    })
  }

  #[expect(clippy::integer_division, reason = "integer division is the intended floor division")]
  pub fn negate_arg_inplace(&mut self) -> Result<(), Report> {
    self.grid = self.grid.negated()?;
    let n = self.y.len();
    for i in 0..n / 2 {
      self.y.swap(i, n - 1 - i);
    }
    Ok(())
  }

  pub fn resample_with_extrap(
    &self,
    grid: &Grid,
    left_extrap: BoundaryBehavior,
    right_extrap: BoundaryBehavior,
  ) -> Result<Self, Report> {
    let n_points = grid.n_points();
    let y_new = (0..n_points)
      .map(|i| self.interp_with_extrap(grid.x_at(i), left_extrap, right_extrap))
      .collect::<Result<Vec<f64>, Report>>()?;
    Self::from_grid_values(*grid, Array1::from_vec(y_new))
  }

  pub fn resample_range_dx_clamped(
    &self,
    x_range: (f64, f64),
    dx: f64,
    anchor: Side,
    max_points: MaxGridPoints,
  ) -> Result<Self, Report> {
    let grid = Grid::from_range_dx_anchored(x_range.0, x_range.1, dx, anchor, max_points)?;
    let x_min = self.grid.x_min();
    let x_max = self.grid.x_max();
    let y_new = (0..grid.n_points())
      .map(|i| {
        self.interp_with_extrap(
          grid.x_at(i).max(x_min).min(x_max),
          BoundaryBehavior::Error,
          BoundaryBehavior::Error,
        )
      })
      .collect::<Result<Vec<f64>, Report>>()?;
    Self::from_grid_values(grid, Array1::from_vec(y_new))
  }
}

struct HardCellEnds {
  boundary: usize,
  anchor: usize,
  side: Side,
}
