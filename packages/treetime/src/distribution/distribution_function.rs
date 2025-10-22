use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use treetime_convolution::{GridFn, InterpElem};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DistributionFunction<T: InterpElem> {
  grid_fn: GridFn<T>,
}

impl<T: InterpElem> DistributionFunction<T> {
  pub fn new(x: Array1<T>, y: Array1<T>) -> Result<Self, Report> {
    assert_eq!(x.shape(), y.shape());
    let grid_fn = GridFn::new(x, y)?;
    Ok(Self { grid_fn })
  }

  pub fn t(&self) -> &Array1<T> {
    self.grid_fn.x()
  }

  pub fn t_mut(&mut self) -> &mut Array1<T> {
    self.grid_fn.x_mut()
  }

  pub fn y(&self) -> &Array1<T> {
    self.grid_fn.y()
  }

  pub fn y_mut(&mut self) -> &mut Array1<T> {
    self.grid_fn.y_mut()
  }

  pub fn interp(&self, x: T) -> Result<T, Report> {
    self.grid_fn.interp(x)
  }

  pub fn interp_many(&self, xs: &Array1<T>) -> Result<Array1<T>, Report> {
    self.grid_fn.interp_many(xs)
  }

  /// Find the most likely time point (x-value corresponding to minimum y-value)
  /// Returns the x-coordinate where the function reaches its minimum y-value
  pub fn likely_time(&self) -> Option<T> {
    let y_values = self.y();
    if let Ok(min_idx) = y_values.argmin() {
      Some(self.t()[min_idx])
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
    let x = array![0.0, 1.0, 2.0, 3.0];
    let y = array![0.0, 10.0, 20.0, 30.0];
    let f = DistributionFunction::new(x, y)?;
    pretty_assert_ulps_eq!(f.interp(1.5)?, 15.0);
    Ok(())
  }

  #[test]
  fn test_distribution_function_interpolate_many() -> Result<(), Report> {
    let x = array![0.0, 1.0, 2.0, 3.0];
    let y = array![0.0, 10.0, 20.0, 30.0];

    let f = DistributionFunction::new(x, y)?;

    let query = array![1.5, 2.25, 2.7];
    let expected = array![15.0, 22.5, 27.0];

    pretty_assert_ulps_eq!(f.interp_many(&query)?, expected);
    Ok(())
  }
}
