use crate::InterpElem;
use crate::grid_iter::GridIter;
use crate::hard_approach_law::Side;
use crate::max_grid_points::MaxGridPoints;
use eyre::Report;
use ndarray::Array1;
use num_traits::Float;
use std::cmp::min;
use std::ops::RangeInclusive;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::float_compare::within_epsilons_of_scale;
use treetime_utils::make_error;

const MAX_COINCIDENCE_ULPS: f64 = 4.0;

const MIN_SPACING_EPSILONS: f64 = 4.0;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Grid<T: InterpElem> {
  x_min: T,
  x_max: T,
  dx: T,
  n_points: usize,
}

impl<T: InterpElem> Grid<T> {
  pub fn from_range_n_points(x_min: T, x_max: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    if x_max <= x_min {
      return make_error!("x_max ({x_max:?}) must be greater than x_min ({x_min:?})");
    }
    let dx = (x_max - x_min) / from_count::<T>(n_points - 1);
    let resolution = <T as From<f64>>::from(MIN_SPACING_EPSILONS) * T::epsilon() * x_min.abs().max(x_max.abs());
    if dx.is_nan() || dx <= resolution {
      return make_error!(
        "Grid of {n_points} points over [{x_min:?}, {x_max:?}] has spacing {dx:?}, not above the rounding resolution {resolution:?} of its coordinates"
      );
    }
    Ok(Self {
      x_min,
      x_max,
      dx,
      n_points,
    })
  }

  pub(crate) fn from_start_dx(x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    if dx <= T::zero() {
      return make_error!("Grid spacing dx ({dx:?}) must be positive");
    }
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    Self::from_range_n_points(x_min, x_min + dx * from_count::<T>(n_points - 1), n_points)
  }

  pub fn from_range_dx(x_min: T, x_max: T, dx: T, max_points: MaxGridPoints) -> Result<Self, Report>
  where
    T: Float,
  {
    let n_points = Self::range_dx_point_count(x_min, x_max, dx, max_points)?;
    Self::from_start_dx(x_min, dx, n_points)
  }

  pub(crate) fn from_range_dx_anchored(
    x_min: T,
    x_max: T,
    dx: T,
    anchor: Side,
    max_points: MaxGridPoints,
  ) -> Result<Self, Report>
  where
    T: Float,
  {
    let n_points = Self::range_dx_point_count(x_min, x_max, dx, max_points)?;
    match anchor {
      Side::Left => Self::from_start_dx(x_min, dx, n_points),
      Side::Right => Self::from_range_n_points(x_max - dx * from_count::<T>(n_points - 1), x_max, n_points),
    }
  }

  pub(crate) fn from_array(x: &Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let n_points = x.len();
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    let grid = Self::from_range_n_points(x[0], x[n_points - 1], n_points)?;
    let scale: f64 = x[0].abs().max(x[n_points - 1].abs()).into();
    match x
      .iter()
      .zip(grid.iter())
      .position(|(&given, at)| !coordinates_coincide(given.into(), at.into(), scale))
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

  fn range_dx_point_count(x_min: T, x_max: T, dx: T, max_points: MaxGridPoints) -> Result<usize, Report>
  where
    T: Float,
  {
    if dx <= T::zero() {
      return make_error!("Grid spacing dx ({dx:?}) must be positive");
    }
    if x_max <= x_min {
      return make_error!("x_max ({x_max:?}) must be greater than x_min ({x_min:?})");
    }
    let required = ((x_max - x_min) / dx + T::one()).round();
    let n_points = max_points.point_count(required.into(), (x_min.into(), x_max.into()), dx.into())?;
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    Ok(n_points)
  }

  pub(crate) fn shifted(&self, delta: T) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::from_range_n_points(self.x_min + delta, self.x_max + delta, self.n_points)
  }

  pub(crate) fn negated(&self) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::from_range_n_points(T::zero() - self.x_max, T::zero() - self.x_min, self.n_points)
  }

  pub(crate) fn scaled(&self, factor: T) -> Result<Self, Report>
  where
    T: Float,
  {
    Self::from_range_n_points(self.x_min * factor, self.x_max * factor, self.n_points)
  }

  pub fn slice(&self, points: RangeInclusive<usize>) -> Result<Self, Report>
  where
    T: Float,
  {
    let (first, last) = (*points.start(), *points.end());
    if last >= self.n_points || first >= last {
      return make_error!(
        "Points {first}..={last} are not a grid of at least 2 points within the grid of {} points",
        self.n_points
      );
    }
    Self::from_range_n_points(self.x_at(first), self.x_at(last), last - first + 1)
  }

  pub fn x_min(&self) -> T {
    self.x_min
  }

  pub fn x_max(&self) -> T {
    self.x_max
  }

  pub fn dx(&self) -> T {
    self.dx
  }

  pub fn coincides(&self, a: f64, b: f64) -> bool
  where
    T: Float,
  {
    let x_min: f64 = self.x_min.into();
    let x_max: f64 = self.x_max.into();
    coordinates_coincide(a, b, x_min.abs().max(x_max.abs()))
  }

  pub fn n_points(&self) -> usize {
    self.n_points
  }

  pub(crate) fn len(&self) -> usize {
    self.n_points
  }

  pub fn x_at(&self, idx: usize) -> T
  where
    T: Float,
  {
    if idx + 1 == self.n_points {
      self.x_max
    } else {
      self.x_min + self.dx * from_count::<T>(idx)
    }
  }

  pub(crate) fn find_interval_index(&self, x: T) -> usize
  where
    T: Float,
  {
    let rel = (x - self.x_min) / self.dx;
    if rel < T::zero() {
      return 0;
    }
    let last_interval = self.n_points - 2;
    rel
      .floor()
      .to_usize()
      .map_or(last_interval, |idx| min(idx, last_interval))
  }

  pub fn to_array(&self) -> Array1<T>
  where
    T: Float,
  {
    self.iter().collect()
  }

  pub fn iter(&self) -> GridIter<T>
  where
    T: Float,
  {
    GridIter::new(*self)
  }
}

pub fn coordinates_coincide(a: f64, b: f64, scale: f64) -> bool {
  within_epsilons_of_scale(a, b, scale, MAX_COINCIDENCE_ULPS)
}

pub(crate) fn from_count<T: InterpElem>(count: usize) -> T {
  count.to_f64_exact().into()
}
