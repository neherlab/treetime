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
  fn test_product_multiplies_the_function_group_with_the_other_factors() {
    let a = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let range = Distribution::range((1.0, 3.0), 0.5);
    let b = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_product(&[&a, &range, &b], MaxGridPoints::default()).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    pretty_assert_ulps_eq!(array![1.0, 2.0, 3.0], actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(
      array![4.5, 5.5, 6.5],
      actual.y(),
      max_ulps = 4,
      "the function costs 1 + t and 2 add to 3 + t, and the range adds its cost 0.5 on its support [1, 3]"
    );
  }

  #[rustfmt::skip]
  #[rstest]
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
  #[trace]
  fn test_product_degenerate_overlap_collapses_to_point_or_empty(
    #[case] (factors, expected): (Vec<Distribution>, Distribution),
  ) {
    let factors = factors.iter().collect_vec();
    let actual = distribution_product(&factors, MaxGridPoints::default()).unwrap();
    assert_eq!(expected, actual);
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
