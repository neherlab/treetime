use crate::make_error;
use eyre::Report;
use num::{FromPrimitive, Zero};
use std::any::type_name;
use std::ops::{Add, Div};

pub trait MeanByKey<T>: Iterator<Item = T> + Sized {
  fn mean_by_key<F, U>(self, mut f: F) -> Result<U, Report>
  where
    F: FnMut(T) -> U,
    U: Add<Output = U> + Div<U, Output = U> + Zero + FromPrimitive + Copy,
  {
    let (sum, count) = self.fold((U::zero(), 0_usize), |(sum, count), item| (sum + f(item), count + 1));
    if count == 0 {
      return make_error!("Cannot compute the mean of an empty sequence");
    }
    let Some(divisor) = U::from_usize(count) else {
      return make_error!(
        "Cannot compute the mean of {count} values: the count does not fit in the key type {}",
        type_name::<U>()
      );
    };
    Ok(sum / divisor)
  }
}

impl<I, T> MeanByKey<T> for I where I: Iterator<Item = T> {}
