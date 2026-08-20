#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_ops::multiply::distribution_multiplication;
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

  /// The two-operand product equals the pairwise multiplication, since both route through the one
  /// co-location primitive. This is the regression guard for the divergence between the endpoint-
  /// inclusive `linspace` grid (pairwise) and the old `resample_range_dx` grid (N-ary): with different
  /// operand spacings the grids used to disagree, so `a * b` differed from `product([a, b])`.
  #[test]
  fn test_product_two_factors_equals_pairwise_multiplication() {
    // Different spacings (dx = 1.0 and dx = 0.5) so the co-located grid is a non-trivial resample.
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();

    let pairwise = distribution_multiplication(&a, &b).unwrap();
    let nary = distribution_product(&[&a, &b]).unwrap();
    // Bit-identical by construction: both call the same `multiply_functions` primitive.
    assert_eq!(pairwise, nary);
  }

  /// Two functions whose supports meet at a single point collapse to a point mass, not an internal
  /// error. The old N-ary fold hit `left >= right` at endpoint contact and raised an internal error;
  /// the pairwise arm already produced a point, and both now share that trichotomy.
  #[test]
  fn test_product_endpoint_contact_returns_point() {
    let a = Distribution::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let b = Distribution::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_product(&[&a, &b]).unwrap();
    // Oracle: v0 converts the one surviving knot at t=1 to a delta; amplitude a(1) * b(1) = 3 * 5.
    let expected = Distribution::point(1.0, 15.0);
    assert_eq!(expected, actual);
  }

  /// Three functions whose hard domains intersect at a single point collapse to a point mass. Exercises
  /// the N-ary single-point case with more than two operands.
  #[test]
  fn test_product_three_functions_single_point_intersection_returns_point() {
    let a = Distribution::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap(); // [0, 1]
    let b = Distribution::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap(); // [1, 2]
    let c = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 1.0]).unwrap(); // [0, 2]

    // Hard intersection: lo = max(0, 1, 0) = 1, hi = min(1, 2, 2) = 1, so the supports touch at t = 1.
    let actual = distribution_product(&[&a, &b, &c]).unwrap();
    // Oracle: amplitude a(1) * b(1) * c(1) = 3 * 5 * 4 (order-independent).
    let expected = Distribution::point(1.0, 60.0);
    assert_eq!(expected, actual);
  }

  /// Functions with disjoint hard (default `Error`) domains produce a legitimate `Empty`: only hard
  /// sides can separate domains, and here they do not overlap.
  #[test]
  fn test_product_disjoint_hard_functions_returns_empty() {
    let a = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let b = Distribution::function(array![5.0, 6.0, 7.0], array![3.0, 2.0, 1.0]).unwrap();

    assert_eq!(Distribution::Empty, distribution_product(&[&a, &b]).unwrap());
  }

  /// Hard support bounds land exactly on grid nodes, because the co-located grid is built with
  /// `Array1::linspace` over the exact intersection. This is what keeps a mode sitting on a hard edge
  /// representable (the endpoint-inclusive grid the recorded contract mandates).
  #[test]
  fn test_product_hard_bounds_land_on_grid_nodes() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![1.0, 1.5, 2.0, 2.5, 3.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    // Intersection [1, 3], both hard; the first and last grid nodes are exactly the hard bounds.
    let Distribution::Function(f) = distribution_product(&[&a, &b]).unwrap() else {
      panic!("Expected Function variant");
    };
    let t = f.t();
    pretty_assert_ulps_eq!(1.0, t[0], max_ulps = 4);
    pretty_assert_ulps_eq!(3.0, t[t.len() - 1], max_ulps = 4);
  }

  /// With operands on different grids, every permutation folds to a bit-identical result. The canonical
  /// operand order removes the floating-point non-associativity of the ordinate sum, so order
  /// independence is exact rather than only within rounding tolerance.
  #[test]
  fn test_product_bit_identical_across_operand_order() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();
    let c = Distribution::function(array![0.5, 1.25, 2.0, 2.75, 3.5], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();

    let abc = distribution_product(&[&a, &b, &c]).unwrap();
    let cba = distribution_product(&[&c, &b, &a]).unwrap();
    let bca = distribution_product(&[&b, &c, &a]).unwrap();
    assert_eq!(abc, cba);
    assert_eq!(abc, bca);
  }
}
