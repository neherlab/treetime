#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_ops::product::distribution_product;
  use itertools::Itertools;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_grid::MaxGridPoints;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_product_three_functions_adds_the_costs_on_a_shared_grid() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();
    let c = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_product(&[&a, &b, &c], MaxGridPoints::default()).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    let expected_y = array![
      1.0 + 5.0 + 2.0,
      2.0 + 4.0 + 2.0,
      3.0 + 3.0 + 2.0,
      4.0 + 2.0 + 2.0,
      5.0 + 1.0 + 2.0
    ];
    pretty_assert_ulps_eq!(array![0.0, 1.0, 2.0, 3.0, 4.0], actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_product_independent_of_factor_order() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();
    let c = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 3.0, 5.0, 3.0, 2.0]).unwrap();

    let Distribution::Function(abc) = distribution_product(&[&a, &b, &c], MaxGridPoints::default()).unwrap() else {
      panic!("Expected Function variant");
    };
    let Distribution::Function(cab) = distribution_product(&[&c, &a, &b], MaxGridPoints::default()).unwrap() else {
      panic!("Expected Function variant");
    };
    pretty_assert_ulps_eq!(abc.y(), cab.y(), max_ulps = 8);
  }

  #[test]
  fn test_product_point_factor_samples_function() {
    let function = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();
    let point = Distribution::point(1.0, 2.0);

    let actual = distribution_product(&[&function, &point], MaxGridPoints::default()).unwrap();
    let expected = Distribution::point(1.0, 2.0 + 2.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_product_two_factors_add_their_costs_on_the_overlap_at_the_finer_spacing() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();

    let actual = distribution_product(&[&a, &b], MaxGridPoints::default()).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    pretty_assert_ulps_eq!(array![1.2, 1.7, 2.2, 2.7], actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(
      array![4.2, 5.7, 7.2, 8.7],
      actual.y(),
      max_ulps = 4,
      "the cost 1 + t of the first factor plus the cost 2 + 2 (t - 1.2) of the second, which is 3 t + 0.6 on the overlap [1.2, 2.7]"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::endpoint_contact_of_two((
    vec![
      Distribution::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap(),
      Distribution::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap(),
    ],
    Distribution::point(1.0, 3.0 + 5.0),
  ))]
  #[case::single_point_intersection_of_three((
    vec![
      Distribution::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap(),
      Distribution::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap(),
      Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 1.0]).unwrap(),
    ],
    Distribution::point(1.0, 3.0 + 5.0 + 4.0),
  ))]
  #[case::empty_factor((
    vec![
      Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
      Distribution::Empty,
      Distribution::function(array![0.0, 1.0, 2.0], array![3.0, 2.0, 1.0]).unwrap(),
    ],
    Distribution::Empty,
  ))]
  #[case::disjoint_hard_functions((
    vec![
      Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
      Distribution::function(array![5.0, 6.0, 7.0], array![3.0, 2.0, 1.0]).unwrap(),
    ],
    Distribution::Empty,
  ))]
  #[trace]
  fn test_product_degenerate_overlap_collapses_to_point_or_empty(
    #[case] (factors, expected): (Vec<Distribution>, Distribution),
  ) {
    let factors = factors.iter().collect_vec();
    let actual = distribution_product(&factors, MaxGridPoints::default()).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_product_hard_bounds_land_on_grid_nodes() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![1.0, 1.5, 2.0, 2.5, 3.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let Distribution::Function(f) = distribution_product(&[&a, &b], MaxGridPoints::default()).unwrap() else {
      panic!("Expected Function variant");
    };
    let t = f.t();
    pretty_assert_ulps_eq!(1.0, t[0], max_ulps = 4);
    pretty_assert_ulps_eq!(3.0, t[t.len() - 1], max_ulps = 4);
  }

  #[test]
  fn test_product_bit_identical_across_operand_order() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let b = Distribution::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();
    let c = Distribution::function(array![0.5, 1.25, 2.0, 2.75, 3.5], array![5.0, 4.0, 3.0, 2.0, 1.0]).unwrap();

    let abc = distribution_product(&[&a, &b, &c], MaxGridPoints::default()).unwrap();
    let cba = distribution_product(&[&c, &b, &a], MaxGridPoints::default()).unwrap();
    let bca = distribution_product(&[&b, &c, &a], MaxGridPoints::default()).unwrap();
    assert_eq!(abc, cba);
    assert_eq!(abc, bca);
  }
}
