#[cfg(test)]
mod __tests__;

use num::{FromPrimitive, Zero};
use std::ops::{Add, Div};

pub trait MeanByKey<T>: Iterator<Item = T> + Sized {
  /// Calculate the arithmetic mean of values extracted from an iterator
  fn mean_by_key<F, U>(self, mut f: F) -> U
  where
    F: FnMut(T) -> U,
    U: Add<Output = U> + Div<U, Output = U> + Zero + FromPrimitive + Copy,
  {
    let (sum, count) = self.fold((U::zero(), 0_usize), |(sum, count), item| (sum + f(item), count + 1));
    if count == 0 {
      U::zero()
    } else {
      sum / U::from_usize(count).unwrap_or_else(U::zero)
    }
  }
}

impl<I, T> MeanByKey<T> for I where I: Iterator<Item = T> {}
