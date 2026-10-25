use crate::InterpElem;
use crate::grid_iter::GridIter;
use crate::hard_approach_law::Side;
use crate::max_grid_points::MaxGridPoints;
use eyre::Report;
use ndarray::Array1;
use num_traits::Float;
use std::cmp::min;
use treetime_utils::array::ndarray::has_uniform_spacing;
use treetime_utils::cast::ToF64Exact;
use treetime_utils::float_compare::within_epsilons_of_scale;
use treetime_utils::make_error;

const MAX_COINCIDENCE_ULPS: f64 = 4.0;

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
  dx: T,
  n_points: usize,
}

impl<T: InterpElem> Grid<T> {
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
    Ok(Self { x_min, dx, n_points })
  }

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
    Ok(Self { x_min, dx, n_points })
  }

  pub fn from_range_dx(x_min: T, x_max: T, dx: T, max_points: MaxGridPoints) -> Result<Self, Report>
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
    Ok(Self { x_min, dx, n_points })
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
    let grid = Self::from_range_dx(x_min, x_max, dx, max_points)?;
    Ok(match anchor {
      Side::Left => grid,
      Side::Right => Self {
        x_min: x_max - dx * from_count::<T>(grid.n_points - 1),
        ..grid
      },
    })
  }

  pub(crate) fn from_array(x: &Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    if x.len() < 2 {
      return make_error!("Grid must have at least 2 points, got {}", x.len());
    }
    if !has_uniform_spacing(x) {
      return make_error!("x array must be uniformly spaced");
    }
    let dx = x[1] - x[0];
    Ok(Self {
      x_min: x[0],
      dx,
      n_points: x.len(),
    })
  }

  pub(crate) fn x_min(&self) -> T {
    self.x_min
  }

  pub(crate) fn x_max(&self) -> T
  where
    T: Float,
  {
    self.x_min + self.dx * from_count::<T>(self.n_points - 1)
  }

  pub(crate) fn dx(&self) -> T {
    self.dx
  }

  pub fn coincides(&self, a: f64, b: f64) -> bool
  where
    T: Float,
  {
    let x_min: f64 = self.x_min.into();
    let x_max: f64 = self.x_max().into();
    coordinates_coincide(a, b, x_min.abs().max(x_max.abs()))
  }

  pub fn n_points(&self) -> usize {
    self.n_points
  }

  pub(crate) fn len(&self) -> usize {
    self.n_points
  }

  pub(crate) fn x_at(&self, idx: usize) -> T
  where
    T: Float,
  {
    self.x_min + self.dx * from_count::<T>(idx)
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

  fn iter(&self) -> GridIter<T>
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
