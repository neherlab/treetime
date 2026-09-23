#[cfg(test)]
mod tests {
  use crate::DistributionPlain;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::distribution_ops::product::distribution_product;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_product_three_functions_is_elementwise_product_on_shared_grid() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();
    let c = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_product(&[&a, &b, &c]).unwrap();
    let DistributionPlain::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    let expected_y = array![10.0, 16.0, 18.0, 16.0, 10.0];
    pretty_assert_ulps_eq!(array![0.0, 1.0, 2.0, 3.0, 4.0], actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_product_independent_of_factor_order() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();
    let c = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 3.0, 5.0, 3.0, 2.0]).unwrap();

    let DistributionPlain::Function(abc) = distribution_product(&[&a, &b, &c]).unwrap() else {
      panic!("Expected Function variant");
    };
    let DistributionPlain::Function(cab) = distribution_product(&[&c, &a, &b]).unwrap() else {
      panic!("Expected Function variant");
    };
    pretty_assert_ulps_eq!(abc.y(), cab.y(), max_ulps = 8);
  }

  #[test]
  fn test_product_empty_factor_returns_empty() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let b = DistributionPlain::function(array![0.0, 1.0, 2.0], array![3.0, 2.0, 1.0]).unwrap();

    let actual = distribution_product(&[&a, &DistributionPlain::Empty, &b]).unwrap();
    assert_eq!(DistributionPlain::Empty, actual);
  }

  #[test]
  fn test_product_point_factor_samples_function() {
    let function = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let point = DistributionPlain::point(1.0, 2.0);

    let actual = distribution_product(&[&function, &point]).unwrap();
    let expected = DistributionPlain::point(1.0, 4.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_product_two_factors_equals_pairwise_multiplication() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = DistributionPlain::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();

    let pairwise = distribution_multiplication(&a, &b).unwrap();
    let nary = distribution_product(&[&a, &b]).unwrap();
    assert_eq!(pairwise, nary);
  }

  #[test]
  fn test_product_endpoint_contact_returns_point() {
    let a = DistributionPlain::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let b = DistributionPlain::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_product(&[&a, &b]).unwrap();
    let expected = DistributionPlain::point(1.0, 15.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_product_three_functions_single_point_intersection_returns_point() {
    let a = DistributionPlain::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let b = DistributionPlain::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();
    let c = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 1.0]).unwrap();

    let actual = distribution_product(&[&a, &b, &c]).unwrap();
    let expected = DistributionPlain::point(1.0, 60.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_product_disjoint_hard_functions_returns_empty() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let b = DistributionPlain::function(array![5.0, 6.0, 7.0], array![3.0, 2.0, 1.0]).unwrap();

    assert_eq!(DistributionPlain::Empty, distribution_product(&[&a, &b]).unwrap());
  }

  #[test]
  fn test_product_hard_bounds_land_on_grid_nodes() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = DistributionPlain::function(array![1.0, 1.5, 2.0, 2.5, 3.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let DistributionPlain::Function(f) = distribution_product(&[&a, &b]).unwrap() else {
      panic!("Expected Function variant");
    };
    let t = f.t();
    pretty_assert_ulps_eq!(1.0, t[0], max_ulps = 4);
    pretty_assert_ulps_eq!(3.0, t[t.len() - 1], max_ulps = 4);
  }

  #[test]
  fn test_product_bit_identical_across_operand_order() {
    let a = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = DistributionPlain::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();
    let c = DistributionPlain::function(array![0.5, 1.25, 2.0, 2.75, 3.5], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();

    let abc = distribution_product(&[&a, &b, &c]).unwrap();
    let cba = distribution_product(&[&c, &b, &a]).unwrap();
    let bca = distribution_product(&[&b, &c, &a]).unwrap();
    assert_eq!(abc, cba);
    assert_eq!(abc, bca);
  }
}
