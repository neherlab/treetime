use approx::UlpsEq;
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use num::Float;
use serde::{Deserialize, Serialize};
use treetime_convolution::grid::Grid;
use treetime_convolution::{GridFn, InterpElem};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DistributionFunction<T: InterpElem> {
  grid_fn: GridFn<T>,
}

impl<T: InterpElem> DistributionFunction<T> {
  pub fn from_arrays(x: &Array1<T>, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays(x, y)?;
    Ok(Self { grid_fn })
  }

  pub fn from_arrays_nonuniform(x: &Array1<T>, y: &Array1<T>) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = GridFn::from_arrays_nonuniform(x, y)?;
    Ok(Self { grid_fn })
  }

  pub fn from_range_values(x_range: (T, T), y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::from_range_values(x_range, y)?;
    Ok(Self { grid_fn })
  }

  pub fn from_start_dx_values(x_min: T, dx: T, y: Array1<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::from_start_dx_values(x_min, dx, y)?;
    Ok(Self { grid_fn })
  }

  pub fn from_n_points<F>((x_min, x_max): (T, T), n_points: usize, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid_fn = GridFn::from_n_points((x_min, x_max), n_points, y_fn)?;
    Ok(Self { grid_fn })
  }

  pub fn from_grid<F>((x_min, x_max): (T, T), dx: T, y_fn: F) -> Result<Self, Report>
  where
    T: Float,
    F: Fn(T) -> T,
  {
    let grid_fn = GridFn::from_grid((x_min, x_max), dx, y_fn)?;
    Ok(Self { grid_fn })
  }

  pub fn constant((x_min, x_max): (T, T), n_points: usize, value: T) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::constant((x_min, x_max), n_points, value)?;
    Ok(Self { grid_fn })
  }

  pub fn zeros((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::zeros((x_min, x_max), n_points)?;
    Ok(Self { grid_fn })
  }

  pub fn ones((x_min, x_max): (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = GridFn::ones((x_min, x_max), n_points)?;
    Ok(Self { grid_fn })
  }

  pub fn t(&self) -> Array1<T>
  where
    T: Float,
  {
    self.grid_fn.x()
  }

  pub fn x_min(&self) -> T {
    self.grid_fn.x_min()
  }

  pub fn x_max(&self) -> T
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

  pub fn interp(&self, x: T) -> Result<T, Report>
  where
    T: Float,
  {
    self.grid_fn.interp(x)
  }

  pub fn interp_many(&self, xs: &Array1<T>) -> Result<Array1<T>, Report>
  where
    T: Float,
  {
    self.grid_fn.interp_many(xs)
  }

  pub fn resample(&self, grid: &Grid<T>) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = self.grid_fn.resample(grid)?;
    Ok(Self { grid_fn })
  }

  pub fn resample_start_dx(&self, x_min: T, dx: T, n_points: usize) -> Result<Self, Report>
  where
    T: Float,
  {
    let grid_fn = self.grid_fn.resample_start_dx(x_min, dx, n_points)?;
    Ok(Self { grid_fn })
  }

  pub fn resample_range_n_points(&self, x_range: (T, T), n_points: usize) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_range_n_points(x_range, n_points)?;
    Ok(Self { grid_fn })
  }

  pub fn resample_range_dx(&self, x_range: (T, T), dx: T) -> Result<Self, Report>
  where
    T: Float + UlpsEq,
  {
    let grid_fn = self.grid_fn.resample_range_dx(x_range, dx)?;
    Ok(Self { grid_fn })
  }

  pub fn len(&self) -> usize {
    self.grid_fn.len()
  }

  pub fn is_empty(&self) -> bool {
    self.grid_fn.len() == 0
  }

  #[must_use]
  pub fn negate_arg(&self) -> Self
  where
    T: Float,
  {
    let grid_fn = self.grid_fn.negate_arg();
    Self { grid_fn }
  }

  pub fn negate_arg_inplace(&mut self)
  where
    T: Float,
  {
    self.grid_fn.negate_arg_inplace();
  }

  /// Find the most likely time point (x-value corresponding to maximum y-value)
  /// Returns the x-coordinate where the function reaches its maximum y-value
  pub fn likely_time(&self) -> Option<T>
  where
    T: Float,
  {
    let y_values = self.y();
    if let Ok(max_idx) = y_values.argmax() {
      Some(self.t()[max_idx])
    } else {
      None
    }
  }
}

#[cfg(test)]
mod tests {
  use crate::distribution::distribution_function::DistributionFunction;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;

  #[test]
  fn test_distribution_function_interpolate() -> Result<(), Report> {
    let f = DistributionFunction::from_range_values((0.0, 3.0), array![0.0, 10.0, 20.0, 30.0])?;
    pretty_assert_ulps_eq!(f.interp(1.5)?, 15.0);
    Ok(())
  }

  #[test]
  fn test_distribution_function_interpolate_many() -> Result<(), Report> {
    let f = DistributionFunction::from_range_values((0.0, 3.0), array![0.0, 10.0, 20.0, 30.0])?;

    let query = array![1.5, 2.25, 2.7];
    let expected = array![15.0, 22.5, 27.0];

    pretty_assert_ulps_eq!(f.interp_many(&query)?, expected);
    Ok(())
  }
}
