use crate::grid_iter::GridIter;
use crate::hard_approach_law::Side;
use crate::max_grid_points::MaxGridPoints;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use num_traits::ToPrimitive;
use std::cmp::min;
use std::fmt::Arguments;
use std::ops::RangeInclusive;
use treetime_utils::cast::{F64_EXACT_INTEGER_LIMIT, ToF64Exact, usize_to_u64};
use treetime_utils::float_compare::within_epsilons_of_scale;
use treetime_utils::make_error;

const MAX_COINCIDENCE_ULPS: f64 = 4.0;

const MAX_ARRAY_DEVIATION_CELLS: f64 = 0.25;

const MIN_UNCHECKED_SPACING_ROUNDINGS: f64 = 8.0;

impl IntoIterator for &Grid {
  type Item = f64;
  type IntoIter = GridIter;

  fn into_iter(self) -> Self::IntoIter {
    self.iter()
  }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grid {
  x_min: f64,
  x_max: f64,
  dx: f64,
  n_points: usize,
}

impl Grid {
  pub fn from_range_n_points(x_min: f64, x_max: f64, n_points: usize) -> Result<Self, Report> {
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    ensure_finite("x_min", x_min)?;
    ensure_finite("x_max", x_max)?;
    if x_max <= x_min {
      return make_error!("x_max ({x_max:?}) must be greater than x_min ({x_min:?})");
    }
    let span = exact_span(n_points, format_args!("over [{x_min:?}, {x_max:?}]"))?;
    let dx = (x_max - x_min) / span;
    if !dx.is_finite() || dx <= 0.0 {
      return make_error!(
        "Grid of {n_points} points over [{x_min:?}, {x_max:?}] has spacing {dx:?}, which is not a finite positive number"
      );
    }
    let grid = Self {
      x_min,
      x_max,
      dx,
      n_points,
    };
    if grid.spacing_outweighs_rounding() {
      return Ok(grid);
    }
    match grid.iter().tuple_windows().position(|(left, right)| right <= left) {
      Some(i) => make_error!(
        "Grid of {n_points} points over [{x_min:?}, {x_max:?}] with spacing {dx:?} does not increase strictly: x[{}] = {:?} is not above x[{i}] = {:?}",
        i + 1,
        grid.x_at(i + 1),
        grid.x_at(i)
      ),
      None => Ok(grid),
    }
  }

  pub(crate) fn from_start_dx(x_min: f64, dx: f64, n_points: usize) -> Result<Self, Report> {
    ensure_finite("x_min", x_min)?;
    ensure_spacing(dx)?;
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    let span = exact_span(n_points, format_args!("from {x_min:?} with spacing {dx:?}"))?;
    Self::from_range_n_points(x_min, x_min + dx * span, n_points)
  }

  pub(crate) fn from_range_dx_anchored(
    x_min: f64,
    x_max: f64,
    dx: f64,
    anchor: Side,
    max_points: MaxGridPoints,
  ) -> Result<Self, Report> {
    let n_points = Self::range_dx_point_count(x_min, x_max, dx, max_points)?;
    match anchor {
      Side::Left => Self::from_start_dx(x_min, dx, n_points),
      Side::Right => {
        let span = exact_span(n_points, format_args!("over [{x_min:?}, {x_max:?}]"))?;
        Self::from_range_n_points(x_max - dx * span, x_max, n_points)
      },
    }
  }

  pub(crate) fn from_array(x: &Array1<f64>) -> Result<Self, Report> {
    let n_points = x.len();
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    if let Some(i) = x.iter().position(|value| !value.is_finite()) {
      return make_error!("x array is not a grid: x[{i}] is {:?}, which is not finite", x[i]);
    }
    if let Some(i) = x.iter().tuple_windows().position(|(left, right)| right <= left) {
      return make_error!(
        "x array is not a grid: x[{}] = {:?} is not above x[{i}] = {:?}",
        i + 1,
        x[i + 1],
        x[i]
      );
    }
    let grid = Self::from_range_n_points(x[0], x[n_points - 1], n_points)?;
    let scale = x[0].abs().max(x[n_points - 1].abs());
    let max_deviation = grid.dx * MAX_ARRAY_DEVIATION_CELLS;
    match x
      .iter()
      .zip(grid.iter())
      .position(|(&given, at)| !coordinates_coincide(given, at, scale) || (given - at).abs() > max_deviation)
    {
      Some(i) => make_error!(
        "x array is not a uniform grid: x[{i}] is {:?}, but the uniform grid from x[0] = {:?} to x[{}] = {:?} has {:?} there",
        x[i],
        x[0],
        n_points - 1,
        x[n_points - 1],
        grid.x_at(i)
      ),
      None => Ok(grid),
    }
  }

  fn range_dx_point_count(x_min: f64, x_max: f64, dx: f64, max_points: MaxGridPoints) -> Result<usize, Report> {
    ensure_finite("x_min", x_min)?;
    ensure_finite("x_max", x_max)?;
    ensure_spacing(dx)?;
    if x_max <= x_min {
      return make_error!("x_max ({x_max:?}) must be greater than x_min ({x_min:?})");
    }
    let required = ((x_max - x_min) / dx + 1.0).round();
    let n_points = max_points.point_count(required, (x_min, x_max), dx)?;
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    Ok(n_points)
  }

  fn spacing_outweighs_rounding(&self) -> bool {
    let rounding = f64::EPSILON * self.x_min.abs().max(self.x_max.abs()) + f64::MIN_POSITIVE;
    self.dx > MIN_UNCHECKED_SPACING_ROUNDINGS * rounding
  }

  pub(crate) fn shifted(&self, delta: f64) -> Result<Self, Report> {
    Self::from_range_n_points(self.x_min + delta, self.x_max + delta, self.n_points)
  }

  pub(crate) fn negated(&self) -> Result<Self, Report> {
    Self::from_range_n_points(0.0 - self.x_max, 0.0 - self.x_min, self.n_points)
  }

  pub(crate) fn scaled(&self, factor: f64) -> Result<Self, Report> {
    Self::from_range_n_points(self.x_min * factor, self.x_max * factor, self.n_points)
  }

  pub fn slice(&self, points: RangeInclusive<usize>) -> Result<Self, Report> {
    let (first, last) = (*points.start(), *points.end());
    if last >= self.n_points || first >= last {
      return make_error!(
        "Points {first}..={last} are not a grid of at least 2 points within the grid of {} points",
        self.n_points
      );
    }
    Self::from_range_n_points(self.x_at(first), self.x_at(last), last - first + 1)
  }

  pub fn x_min(&self) -> f64 {
    self.x_min
  }

  pub fn x_max(&self) -> f64 {
    self.x_max
  }

  pub fn dx(&self) -> f64 {
    self.dx
  }

  pub fn coincides(&self, a: f64, b: f64) -> bool {
    coordinates_coincide(a, b, self.x_min.abs().max(self.x_max.abs()))
  }

  pub fn n_points(&self) -> usize {
    self.n_points
  }

  pub(crate) fn len(&self) -> usize {
    self.n_points
  }

  pub fn x_at(&self, idx: usize) -> f64 {
    assert!(
      idx < self.n_points,
      "index {idx} is outside the grid of {} points",
      self.n_points
    );
    if idx + 1 == self.n_points {
      self.x_max
    } else {
      self.x_min + self.dx * idx.to_f64_exact()
    }
  }

  pub(crate) fn find_interval_index(&self, x: f64) -> usize {
    let rel = (x - self.x_min) / self.dx;
    if rel < 0.0 {
      return 0;
    }
    let last_interval = self.n_points - 2;
    rel
      .floor()
      .to_usize()
      .map_or(last_interval, |idx| min(idx, last_interval))
  }

  pub fn to_array(&self) -> Array1<f64> {
    self.iter().collect()
  }

  pub fn iter(&self) -> GridIter {
    GridIter::new(*self)
  }
}

pub fn coordinates_coincide(a: f64, b: f64, scale: f64) -> bool {
  within_epsilons_of_scale(a, b, scale, MAX_COINCIDENCE_ULPS)
}

fn ensure_finite(name: &str, value: f64) -> Result<(), Report> {
  if value.is_finite() {
    Ok(())
  } else {
    make_error!("Grid {name} ({value:?}) must be finite")
  }
}

fn ensure_spacing(dx: f64) -> Result<(), Report> {
  if dx.is_finite() && dx > 0.0 {
    Ok(())
  } else {
    make_error!("Grid spacing dx ({dx:?}) must be a finite positive number")
  }
}

fn exact_span(n_points: usize, extent: Arguments<'_>) -> Result<f64, Report> {
  if usize_to_u64(n_points) >= F64_EXACT_INTEGER_LIMIT {
    return make_error!(
      "Grid of {n_points} points {extent} has more points than {}, the largest count with exact coordinates",
      MaxGridPoints::MAX
    );
  }
  Ok((n_points - 1).to_f64_exact())
}
