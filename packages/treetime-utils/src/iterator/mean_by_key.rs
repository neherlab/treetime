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

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_ulps_eq;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_mean_by_key_i32() {
    let data = [1, 2, 3, 4, 5];
    let mean: i32 = data.iter().mean_by_key(|x| *x * 2);
    assert_eq!(mean, 6);
  }

  #[test]
  fn test_mean_by_key_f64() {
    let data = [1.0, 2.0, 3.0, 4.0];
    let mean: f64 = data.iter().mean_by_key(|x| *x);
    assert_ulps_eq!(mean, 2.5);
  }

  #[test]
  fn test_mean_by_key_f32() {
    let data = [1.0_f32, 2.0, 3.0];
    let mean: f32 = data.iter().mean_by_key(|x| *x * 2.0);
    assert_ulps_eq!(mean, 4.0);
  }

  #[test]
  fn test_mean_by_key_empty() {
    let data: Vec<i32> = vec![];
    let mean: i32 = data.iter().mean_by_key(|x| *x);
    assert_eq!(mean, 0);
  }

  #[test]
  fn test_mean_by_key_single() {
    let data = [42.0];
    let mean: f64 = data.iter().mean_by_key(|x| *x);
    assert_ulps_eq!(mean, 42.0);
  }

  #[test]
  fn test_mean_by_key_struct() {
    #[derive(Clone)]
    struct Point {
      x: f64,
      y: f64,
    }

    let points = [
      Point { x: 1.0, y: 2.0 },
      Point { x: 3.0, y: 4.0 },
      Point { x: 5.0, y: 6.0 },
    ];

    let mean_x: f64 = points.iter().mean_by_key(|p| p.x);
    let mean_y: f64 = points.iter().mean_by_key(|p| p.y);

    assert_ulps_eq!(mean_x, 3.0);
    assert_ulps_eq!(mean_y, 4.0);
  }

  #[test]
  fn test_mean_by_key_owned() {
    let data = vec![1, 2, 3, 4];
    let mean: i32 = data.into_iter().mean_by_key(|x| x * 3);
    assert_eq!(mean, 7);
  }

  #[test]
  fn test_mean_by_key_negative_numbers() {
    let data = [-5, -2, -1, 0, 1, 2, 5];
    let mean: i32 = data.iter().mean_by_key(|x| *x);
    assert_eq!(mean, 0);
  }

  #[test]
  fn test_mean_by_key_negative_floats() {
    let data = [-3.5, -1.5, 0.0, 1.5, 3.5];
    let mean: f64 = data.iter().mean_by_key(|x| *x);
    assert_ulps_eq!(mean, 0.0);
  }

  #[test]
  fn test_mean_by_key_large_collection() {
    let data: Vec<i32> = (1..=1000).collect();
    let mean: i32 = data.iter().mean_by_key(|x| *x);
    assert_eq!(mean, 500);
  }

  #[test]
  fn test_mean_by_key_all_zeros() {
    let data = [0, 0, 0, 0, 0];
    let mean: i32 = data.iter().mean_by_key(|x| *x);
    assert_eq!(mean, 0);
  }

  #[test]
  fn test_mean_by_key_u32() {
    let data: Vec<u32> = vec![1, 3, 5, 7, 9];
    let mean: u32 = data.iter().mean_by_key(|x| *x);
    assert_eq!(mean, 5);
  }

  #[test]
  fn test_mean_by_key_i64() {
    let data: Vec<i64> = vec![1000000000, 2000000000, 3000000000];
    let mean: i64 = data.iter().mean_by_key(|x| *x);
    assert_eq!(mean, 2000000000);
  }

  #[test]
  fn test_mean_by_key_complex_transformation() {
    let data = [1, 2, 3, 4];
    let mean: i32 = data.iter().mean_by_key(|x| x * x + 1);
    assert_eq!(mean, 8);
  }

  #[test]
  fn test_mean_by_key_very_small_floats() {
    let data = [1e-10, 2e-10, 3e-10];
    let mean: f64 = data.iter().mean_by_key(|x| *x);
    assert_ulps_eq!(mean, 2e-10);
  }

  #[test]
  fn test_mean_by_key_mixed_precision() {
    let data = [0.1, 0.2, 0.3];
    let mean: f64 = data.iter().mean_by_key(|x| *x);
    assert_ulps_eq!(mean, 0.2);
  }

  #[test]
  fn test_mean_by_key_chained_iterator() {
    let data = [1, 2, 3, 4, 5, 6];
    let mean: i32 = data.iter().filter(|&&x| x % 2 == 0).mean_by_key(|x| *x);
    assert_eq!(mean, 4);
  }

  #[test]
  fn test_mean_by_key_absolute_values() {
    let data = [-3, -1, 1, 3];
    let mean: i32 = data.iter().mean_by_key(|x| if *x < 0 { -*x } else { *x });
    assert_eq!(mean, 2);
  }

  #[test]
  fn test_mean_by_key_option_values() {
    let data = [Some(1), Some(2), None, Some(3)];
    let mean: i32 = data.iter().filter_map(|x| *x).mean_by_key(|x| x);
    assert_eq!(mean, 2);
  }
}
