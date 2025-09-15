use eyre::Report;
use ndarray::{Array1, Ix1, OwnedRepr};
use ndarray_interp::interp1d::{Interp1D, Interp1DBuilder, Linear};
use num_traits::{Num, NumCast};
use std::fmt::Debug;

pub trait InterpElem: Num + NumCast + Debug + Send + PartialOrd + Copy {}

impl InterpElem for f64 {}

/// Represents an arbitrary smooth function
#[derive(Debug)]
pub struct DistributionFunction<T: InterpElem> {
  interp: Interp1D<OwnedRepr<T>, OwnedRepr<T>, Ix1, Linear>,
}

impl<T: InterpElem> PartialEq for DistributionFunction<T> {
  fn eq(&self, other: &Self) -> bool {
    format!("{self:?}") == format!("{other:?}")
  }
}

impl<T: InterpElem> Eq for DistributionFunction<T> {}

impl<T: InterpElem> DistributionFunction<T> {
  pub fn new(x: Array1<T>, y: Array1<T>) -> Result<Self, Report> {
    assert_eq!(x.shape(), y.shape());
    let interp: Interp1D<OwnedRepr<T>, OwnedRepr<T>, Ix1, Linear> = Interp1DBuilder::new(y).x(x).build()?;
    Ok(Self { interp })
  }

  pub fn t(&self) -> &Array1<T> {
    self.interp.x()
  }

  pub fn t_mut(&mut self) -> &mut Array1<T> {
    self.interp.x_mut()
  }

  pub fn y(&self) -> &Array1<T> {
    self.interp.data()
  }

  pub fn y_mut(&mut self) -> &mut Array1<T> {
    self.interp.data_mut()
  }

  pub fn interp(&self, x: T) -> Result<T, Report> {
    let y = self.interp.interp_scalar(x)?;
    Ok(y)
  }

  pub fn interp_many(&self, xs: &Array1<T>) -> Result<Array1<T>, Report> {
    let ys = self.interp.interp_array(xs)?;
    Ok(ys)
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
