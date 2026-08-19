#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_ops::product::distribution_product;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  /// The N-ary product co-locates all `Function` factors on one shared grid and reduces once. On a
  /// common grid (identical `x`, so resampling is the identity) the fold is the exact elementwise
  /// product of the operand ordinates -- an independent oracle for `distribution_product`.
  #[test]
  fn test_product_three_functions_is_elementwise_product_on_shared_grid() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();
    let c = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_product(&[&a, &b, &c]).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    // Oracle: elementwise a * b * c over the shared grid, e.g. 1*5*2, 2*4*2, 3*3*2, 4*2*2, 5*1*2.
    let expected_y = array![10.0, 16.0, 18.0, 16.0, 10.0];
    pretty_assert_ulps_eq!(array![0.0, 1.0, 2.0, 3.0, 4.0], actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  /// The product is a product of independent factors, so its value does not depend on factor order.
  #[test]
  fn test_product_independent_of_factor_order() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();
    let c = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 3.0, 5.0, 3.0, 2.0]).unwrap();

    let Distribution::Function(abc) = distribution_product(&[&a, &b, &c]).unwrap() else {
      panic!("Expected Function variant");
    };
    let Distribution::Function(cab) = distribution_product(&[&c, &a, &b]).unwrap() else {
      panic!("Expected Function variant");
    };
    pretty_assert_ulps_eq!(abc.y(), cab.y(), max_ulps = 8);
  }

  /// An `Empty` factor makes the whole product empty (an empty operand is disjoint from every domain).
  #[test]
  fn test_product_empty_factor_returns_empty() {
    let a = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let b = Distribution::function(array![0.0, 1.0, 2.0], array![3.0, 2.0, 1.0]).unwrap();

    let actual = distribution_product(&[&a, &Distribution::Empty, &b]).unwrap();
    assert_eq!(Distribution::Empty, actual);
  }

  /// A `Point` factor samples the function product at its position, collapsing the result to a point.
  #[test]
  fn test_product_point_factor_samples_function() {
    let function = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let point = Distribution::point(1.0, 2.0);

    let actual = distribution_product(&[&function, &point]).unwrap();
    // Oracle: the point at t=1 samples the function value 2.0 and scales by its amplitude 2.0.
    let expected = Distribution::point(1.0, 4.0);
    assert_eq!(expected, actual);
  }
}
