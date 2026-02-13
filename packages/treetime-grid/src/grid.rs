use crate::InterpElem;
use crate::grid_iter::GridIter;
use approx::UlpsEq;
use eyre::Report;
use ndarray::Array1;
use num_traits::Float;
use serde::{Deserialize, Serialize};
use std::cmp::min;
use treetime_utils::array::ndarray::has_uniform_spacing;
use treetime_utils::make_error;

/// Uniform grid parameters and indexing operations
///
/// Encapsulates grid parameters (x_min, dx, n_points) and provides methods for
/// indexing, coordinate calculation, and iteration.
///
/// # Invariants
///
/// - Grid must be uniformly spaced with spacing `dx`
/// - Grid starts at `x_min`
/// - Must contain at least 2 points
/// - `dx` must be positive
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Grid<T: InterpElem> {
  x_min: T,
  dx: T,
  n_points: usize,
}

impl<T: InterpElem> Grid<T> {
  pub fn from_start_dx(x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
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
    let dx = (x_max - x_min) / T::from(n_points - 1).unwrap();
    Ok(Self { x_min, dx, n_points })
  }

  pub fn from_range_dx(x_min: T, x_max: T, dx: T) -> Result<Self, Report>
  where
    T: Float,
  {
    if dx <= T::zero() {
      return make_error!("Grid spacing dx ({dx:?}) must be positive");
    }
    if x_max <= x_min {
      return make_error!("x_max ({x_max:?}) must be greater than x_min ({x_min:?})");
    }
    let n_points = ((x_max - x_min) / dx + T::one()).round().to_usize().unwrap();
    if n_points < 2 {
      return make_error!("Grid must have at least 2 points, got {n_points}");
    }
    Ok(Self { x_min, dx, n_points })
  }

  pub fn from_array(x: &Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
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

  pub fn x_min(&self) -> T {
    self.x_min
  }

  pub fn x_max(&self) -> T
  where
    T: Float,
  {
    self.x_min + self.dx * T::from(self.n_points - 1).unwrap()
  }

  pub fn x_range(&self) -> (T, T)
  where
    T: Float,
  {
    (self.x_min(), self.x_max())
  }

  pub fn dx(&self) -> T {
    self.dx
  }

  pub fn n_points(&self) -> usize {
    self.n_points
  }

  pub fn len(&self) -> usize {
    self.n_points
  }

  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  /// Computes x coordinate at given index
  pub fn x_at(&self, idx: usize) -> T
  where
    T: Float,
  {
    self.x_min + self.dx * T::from(idx).unwrap()
  }

  /// Finds grid interval index containing given x value
  pub fn find_interval_index(&self, x: T) -> usize
  where
    T: Float,
  {
    let idx = ((x - self.x_min) / self.dx).floor().to_usize().unwrap();
    min(idx, self.n_points - 2)
  }

  /// Generates array of all x coordinates
  pub fn to_array(&self) -> Array1<T>
  where
    T: Float,
  {
    self.iter().collect()
  }

  /// Returns an iterator over grid x coordinates
  pub fn iter(&self) -> GridIter<T>
  where
    T: Float,
  {
    GridIter::new(*self)
  }
}

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
