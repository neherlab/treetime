#[cfg(test)]
mod tests {
  use crate::{DistributionPlain as Distribution, distribution_division, distribution_multiplication};
  use ndarray::{Array1, array};
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_ulps_eq, prop_assert_ulps_eq};

  proptest! {
    #[test]
    fn test_prop_multiply_range_function_preserves_intersection(
      start_hundredths in -200_i16..600,
      width_hundredths in 1_u16..600,
    ) {
      let start = f64::from(start_hundredths) / 100.0;
      let end = start + f64::from(width_hundredths) / 100.0;
      let range = Distribution::range((start, end), 3.0);
      let function = linear_function();

      let actual = distribution_multiplication(&range, &function).unwrap();
      assert_intersection(&actual, (start.max(0.0), end.min(4.0)))?;
    }

    #[test]
    fn test_prop_multiply_range_function_commutative(
      start_hundredths in -200_i16..600,
      width_hundredths in 1_u16..600,
    ) {
      let start = f64::from(start_hundredths) / 100.0;
      let end = start + f64::from(width_hundredths) / 100.0;
      let range = Distribution::range((start, end), 3.0);
      let function = linear_function();

      let range_function = distribution_multiplication(&range, &function).unwrap();
      let function_range = distribution_multiplication(&function, &range).unwrap();
      prop_assert_eq!(range_function, function_range);
    }

    #[test]
    fn test_prop_multiply_function_function_commutative(
      left_points in 2_usize..20,
      right_points in 2_usize..20,
      left_width_hundredths in 100_u16..500,
      right_width_hundredths in 100_u16..500,
    ) {
      let left = linear_function_n_points(left_points, f64::from(left_width_hundredths) / 100.0);
      let right = linear_function_n_points(right_points, f64::from(right_width_hundredths) / 100.0);

      let left_right = distribution_multiplication(&left, &right).unwrap();
      let right_left = distribution_multiplication(&right, &left).unwrap();
      prop_assert_eq!(left_right, right_left);
    }

    #[test]
    fn test_prop_divide_recovers_range_amplitude(
      start_hundredths in 0_u16..400,
      width_hundredths in 1_u16..600,
    ) {
      let start = f64::from(start_hundredths) / 100.0;
      let end = start + f64::from(width_hundredths) / 100.0;
      let range = Distribution::range((start, end), 3.0);
      let function = linear_function();
      let product = distribution_multiplication(&range, &function).unwrap();

      let actual = distribution_division(&product, &function).unwrap();
      let expected = Array1::from_elem(actual.y().len(), 3.0);
      prop_assert_array_ulps_eq!(expected, actual.y(), max_ulps = 4);
      assert_intersection(&actual, (start, end.min(4.0)))?;
    }
  }

  #[allow(clippy::float_cmp)] // Point support requires exact analytical endpoint contact.
  fn assert_intersection(actual: &Distribution, (start, end): (f64, f64)) -> Result<(), TestCaseError> {
    if start > end {
      prop_assert!(matches!(actual, Distribution::Empty));
    } else if start == end {
      let Distribution::Point(point) = actual else {
        return Err(TestCaseError::fail(format!(
          "expected Point at {start}, got {actual:?}"
        )));
      };
      prop_assert_ulps_eq!(start, point.t(), max_ulps = 4);
    } else {
      let Distribution::Function(function) = actual else {
        return Err(TestCaseError::fail(format!(
          "expected Function on [{start}, {end}], got {actual:?}"
        )));
      };
      prop_assert_ulps_eq!(start, function.x_min(), max_ulps = 4);
      prop_assert_ulps_eq!(end, function.x_max(), max_ulps = 4);
    }
    Ok(())
  }

  fn linear_function() -> Distribution {
    Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap()
  }

  fn linear_function_n_points(n_points: usize, width: f64) -> Distribution {
    let x = Array1::linspace(0.0, width, n_points);
    let y = x.mapv(|value| value + 1.0);
    Distribution::function(x, y).unwrap()
  }
}
