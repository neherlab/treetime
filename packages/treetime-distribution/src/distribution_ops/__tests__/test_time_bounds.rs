#[cfg(test)]
mod tests {
  use crate::DistributionPlain;
  use crate::distribution_ops::time_bounds::{
    distribution_support_n_points, distribution_time_bounds_contains, distribution_time_bounds_intersection,
    distribution_time_bounds_overlaps, distribution_time_bounds_union,
  };
  use ndarray::array;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::fractional_ceils_up(     (0.0, 2.4),        1.0, 4)]
  #[case::small_fraction_ceils_up( (0.0, 2.1),        1.0, 4)]
  #[case::exact_multiple_no_extra( (0.0, 3.0),        1.0, 4)]
  #[case::minimum_two(             (0.0, 0.4),        1.0, 2)]
  #[case::maximum_safety_cap(      (0.0, 2_000_000.0), 1.0, 1_000_000)]
  #[trace]
  fn test_distribution_support_n_points_uses_spacing_contract(
    #[case] bounds: (f64, f64),
    #[case] dx: f64,
    #[case] expected: usize,
  ) {
    let actual = distribution_support_n_points(bounds, dx).unwrap();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::two_points_ordered(DistributionPlain::point(1.0, 1.0),           DistributionPlain::point(2.0, 1.0),          (1.0, 2.0)  )]
  #[case::two_points_reversed(DistributionPlain::point(2.0, 1.0),           DistributionPlain::point(1.0, 1.0),          (1.0, 2.0)  )]
  #[case::overlapping_ranges(DistributionPlain::range((1.0, 3.0), 1.0),    DistributionPlain::range((2.0, 4.0), 1.0),   (1.0, 4.0)  )]
  #[case::overlapping_ranges_reversed(DistributionPlain::range((2.0, 4.0), 1.0),    DistributionPlain::range((1.0, 3.0), 1.0),   (1.0, 4.0)  )]
  #[case::disjoint_ranges(DistributionPlain::range((1.0, 2.0), 1.0),    DistributionPlain::range((3.0, 4.0), 1.0),   (1.0, 4.0)  )]
  #[case::nested_ranges(DistributionPlain::range((1.0, 5.0), 1.0),    DistributionPlain::range((2.0, 3.0), 1.0),   (1.0, 5.0)  )]
  #[case::point_and_range(DistributionPlain::point(1.0, 1.0),           DistributionPlain::range((2.0, 3.0), 1.0),   (1.0, 3.0)  )]
  #[case::two_functions(
    DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
    DistributionPlain::function(array![1.5, 2.5, 3.5], array![1.0, 2.0, 3.0]).unwrap(),
    (0.0, 3.5)
  )]
  #[case::same_point(DistributionPlain::point(5.0, 1.0),           DistributionPlain::point(5.0, 1.0),          (5.0, 5.0)  )]
  #[case::negative_and_positive(DistributionPlain::range((-2.0, -1.0), 1.0),  DistributionPlain::range((1.0, 2.0), 1.0),   (-2.0, 2.0) )]
  #[trace]
  fn test_distribution_time_bounds_union(
    #[case] dist_a: DistributionPlain,
    #[case] dist_b: DistributionPlain,
    #[case] expected: (f64, f64),
  ) {
    let actual = distribution_time_bounds_union(&dist_a, &dist_b);
    assert_eq!(Some(expected), actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::overlapping_ranges(DistributionPlain::range((1.0, 3.0), 1.0),   DistributionPlain::range((2.0, 4.0), 1.0),   Some((2.0, 3.0)))]
  #[case::overlapping_ranges_reversed(DistributionPlain::range((2.0, 4.0), 1.0),   DistributionPlain::range((1.0, 3.0), 1.0),   Some((2.0, 3.0)))]
  #[case::disjoint_ranges(DistributionPlain::range((1.0, 2.0), 1.0),   DistributionPlain::range((3.0, 4.0), 1.0),   None            )]
  #[case::disjoint_ranges_reversed(DistributionPlain::range((3.0, 4.0), 1.0),   DistributionPlain::range((1.0, 2.0), 1.0),   None            )]
  #[case::nested_outer_first(DistributionPlain::range((1.0, 5.0), 1.0),   DistributionPlain::range((2.0, 3.0), 1.0),   Some((2.0, 3.0)))]
  #[case::nested_inner_first(DistributionPlain::range((2.0, 3.0), 1.0),   DistributionPlain::range((1.0, 5.0), 1.0),   Some((2.0, 3.0)))]
  #[case::point_inside_range(DistributionPlain::point(2.5, 1.0),          DistributionPlain::range((2.0, 3.0), 1.0),   Some((2.5, 2.5)))]
  #[case::range_contains_point(DistributionPlain::range((2.0, 3.0), 1.0),   DistributionPlain::point(2.5, 1.0),          Some((2.5, 2.5)))]
  #[case::point_outside_range(DistributionPlain::point(1.0, 1.0),          DistributionPlain::range((2.0, 3.0), 1.0),   None            )]
  #[case::two_functions(
    DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
    DistributionPlain::function(array![1.5, 2.5, 3.5], array![1.0, 2.0, 3.0]).unwrap(),
    Some((1.5, 2.0))
  )]
  #[case::adjacent_ranges(DistributionPlain::range((1.0, 2.0), 1.0),   DistributionPlain::range((2.0, 3.0), 1.0),   Some((2.0, 2.0)))]
  #[case::same_point(DistributionPlain::point(5.0, 1.0),          DistributionPlain::point(5.0, 1.0),          Some((5.0, 5.0)))]
  #[case::different_points(DistributionPlain::point(5.0, 1.0),          DistributionPlain::point(6.0, 1.0),          None            )]
  #[trace]
  fn test_distribution_time_bounds_intersection(
    #[case] dist_a: DistributionPlain,
    #[case] dist_b: DistributionPlain,
    #[case] expected: Option<(f64, f64)>,
  ) {
    let actual = distribution_time_bounds_intersection(&dist_a, &dist_b);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::outer_contains_inner(DistributionPlain::range((1.0, 5.0), 1.0),   DistributionPlain::range((2.0, 3.0), 1.0),   true )]
  #[case::inner_not_contains_outer(DistributionPlain::range((2.0, 3.0), 1.0),   DistributionPlain::range((1.0, 5.0), 1.0),   false)]
  #[case::range_contains_point(DistributionPlain::range((1.0, 5.0), 1.0),   DistributionPlain::point(3.0, 1.0),          true )]
  #[case::point_not_contains_range(DistributionPlain::point(3.0, 1.0),          DistributionPlain::range((1.0, 5.0), 1.0),   false)]
  #[case::overlapping_not_contained(DistributionPlain::range((1.0, 3.0), 1.0),   DistributionPlain::range((2.0, 4.0), 1.0),   false)]
  #[case::same_range(DistributionPlain::range((1.0, 5.0), 1.0),   DistributionPlain::range((1.0, 5.0), 1.0),   true )]
  #[case::same_point(DistributionPlain::point(5.0, 1.0),          DistributionPlain::point(5.0, 1.0),          true )]
  #[case::inner_extends_left(DistributionPlain::range((1.0, 5.0), 1.0),   DistributionPlain::range((0.5, 3.0), 1.0),   false)]
  #[case::inner_extends_right(DistributionPlain::range((1.0, 5.0), 1.0),   DistributionPlain::range((3.0, 5.5), 1.0),   false)]
  #[case::function_contains_function(
    DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0], array![1.0, 2.0, 3.0, 2.0, 1.0, 0.5]).unwrap(),
    DistributionPlain::function(array![1.0, 2.0, 3.0], array![1.0, 2.0, 3.0]).unwrap(),
    true
  )]
  #[trace]
  fn test_distribution_time_bounds_contains(
    #[case] outer: DistributionPlain,
    #[case] inner: DistributionPlain,
    #[case] expected: bool,
  ) {
    let actual = distribution_time_bounds_contains(&outer, &inner);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::overlapping_ranges(DistributionPlain::range((1.0, 3.0), 1.0), DistributionPlain::range((2.0, 4.0), 1.0), true )]
  #[case::disjoint_ranges(DistributionPlain::range((1.0, 2.0), 1.0), DistributionPlain::range((3.0, 4.0), 1.0), false)]
  #[case::adjacent_ranges(DistributionPlain::range((1.0, 2.0), 1.0), DistributionPlain::range((2.0, 3.0), 1.0), true )]
  #[case::nested_ranges(DistributionPlain::range((1.0, 5.0), 1.0), DistributionPlain::range((2.0, 3.0), 1.0), true )]
  #[case::point_inside_range(DistributionPlain::point(2.5, 1.0),        DistributionPlain::range((2.0, 3.0), 1.0), true )]
  #[case::point_outside_range(DistributionPlain::point(1.0, 1.0),        DistributionPlain::range((2.0, 3.0), 1.0), false)]
  #[case::same_point(DistributionPlain::point(5.0, 1.0),        DistributionPlain::point(5.0, 1.0),        true )]
  #[case::different_points(DistributionPlain::point(5.0, 1.0),        DistributionPlain::point(6.0, 1.0),        false)]
  #[trace]
  fn test_distribution_time_bounds_overlaps(
    #[case] dist_a: DistributionPlain,
    #[case] dist_b: DistributionPlain,
    #[case] expected: bool,
  ) {
    let actual = distribution_time_bounds_overlaps(&dist_a, &dist_b);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_distribution_time_bounds_union_commutativity() {
    let dist_a = DistributionPlain::range((1.0, 3.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 4.0), 1.0);

    let result_ab = distribution_time_bounds_union(&dist_a, &dist_b);
    let result_ba = distribution_time_bounds_union(&dist_b, &dist_a);

    assert_eq!(result_ab, result_ba);
  }

  #[test]
  fn test_distribution_time_bounds_intersection_commutativity() {
    let dist_a = DistributionPlain::range((1.0, 3.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 4.0), 1.0);

    let result_ab = distribution_time_bounds_intersection(&dist_a, &dist_b);
    let result_ba = distribution_time_bounds_intersection(&dist_b, &dist_a);

    assert_eq!(result_ab, result_ba);
  }

  #[test]
  fn test_distribution_time_bounds_overlaps_symmetry() {
    let dist_a = DistributionPlain::range((1.0, 3.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 4.0), 1.0);

    let result_ab = distribution_time_bounds_overlaps(&dist_a, &dist_b);
    let result_ba = distribution_time_bounds_overlaps(&dist_b, &dist_a);

    assert_eq!(result_ab, result_ba);
  }

  #[test]
  fn test_distribution_time_bounds_union_associativity() {
    let dist_a = DistributionPlain::range((1.0, 2.0), 1.0);
    let dist_b = DistributionPlain::range((3.0, 4.0), 1.0);
    let dist_c = DistributionPlain::range((5.0, 6.0), 1.0);

    let (t_min_ab, t_max_ab) = distribution_time_bounds_union(&dist_a, &dist_b).unwrap();
    let dist_ab = DistributionPlain::range((t_min_ab, t_max_ab), 1.0);
    let result_ab_c = distribution_time_bounds_union(&dist_ab, &dist_c);

    let (t_min_bc, t_max_bc) = distribution_time_bounds_union(&dist_b, &dist_c).unwrap();
    let dist_bc = DistributionPlain::range((t_min_bc, t_max_bc), 1.0);
    let result_a_bc = distribution_time_bounds_union(&dist_a, &dist_bc);

    assert_eq!(result_ab_c, result_a_bc);
  }

  #[test]
  fn test_distribution_time_bounds_intersection_associativity() {
    let dist_a = DistributionPlain::range((1.0, 6.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 5.0), 1.0);
    let dist_c = DistributionPlain::range((3.0, 4.0), 1.0);

    let int_ab = distribution_time_bounds_intersection(&dist_a, &dist_b).unwrap();
    let dist_ab = DistributionPlain::range(int_ab, 1.0);
    let result_ab_c = distribution_time_bounds_intersection(&dist_ab, &dist_c);

    let int_bc = distribution_time_bounds_intersection(&dist_b, &dist_c).unwrap();
    let dist_bc = DistributionPlain::range(int_bc, 1.0);
    let result_a_bc = distribution_time_bounds_intersection(&dist_a, &dist_bc);

    assert_eq!(result_ab_c, result_a_bc);
  }

  #[test]
  fn test_distribution_time_bounds_idempotence() {
    let dist = DistributionPlain::range((1.0, 3.0), 1.0);

    let union_result = distribution_time_bounds_union(&dist, &dist);
    let expected = dist.time_bounds();
    assert_eq!(expected, union_result);

    let intersection_result = distribution_time_bounds_intersection(&dist, &dist).unwrap();
    assert_eq!(expected, Some(intersection_result));
  }

  #[test]
  fn test_distribution_time_bounds_contains_reflexivity() {
    let dist = DistributionPlain::range((1.0, 3.0), 1.0);
    let actual = distribution_time_bounds_contains(&dist, &dist);
    assert!(actual);
  }

  #[test]
  fn test_distribution_time_bounds_contains_transitivity() {
    let dist_outer = DistributionPlain::range((1.0, 10.0), 1.0);
    let dist_middle = DistributionPlain::range((2.0, 8.0), 1.0);
    let dist_inner = DistributionPlain::range((3.0, 7.0), 1.0);

    assert!(distribution_time_bounds_contains(&dist_outer, &dist_middle));
    assert!(distribution_time_bounds_contains(&dist_middle, &dist_inner));
    assert!(distribution_time_bounds_contains(&dist_outer, &dist_inner));
  }

  #[test]
  fn test_distribution_time_bounds_overlaps_reflexivity() {
    let dist = DistributionPlain::range((1.0, 3.0), 1.0);
    let actual = distribution_time_bounds_overlaps(&dist, &dist);
    assert!(actual);
  }

  #[test]
  fn test_distribution_time_bounds_intersection_implies_overlap() {
    let dist_a = DistributionPlain::range((1.0, 3.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 4.0), 1.0);

    let has_intersection = distribution_time_bounds_intersection(&dist_a, &dist_b).is_some();
    let overlaps = distribution_time_bounds_overlaps(&dist_a, &dist_b);

    assert_eq!(has_intersection, overlaps);
  }

  #[test]
  fn test_distribution_time_bounds_contains_implies_overlap() {
    let outer = DistributionPlain::range((1.0, 5.0), 1.0);
    let inner = DistributionPlain::range((2.0, 3.0), 1.0);

    assert!(distribution_time_bounds_contains(&outer, &inner));
    assert!(distribution_time_bounds_overlaps(&outer, &inner));
  }

  #[test]
  fn test_distribution_time_bounds_union_contains_both() {
    let dist_a = DistributionPlain::range((1.0, 3.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 4.0), 1.0);

    let (t_min, t_max) = distribution_time_bounds_union(&dist_a, &dist_b).unwrap();
    let union_dist = DistributionPlain::range((t_min, t_max), 1.0);

    assert!(distribution_time_bounds_contains(&union_dist, &dist_a));
    assert!(distribution_time_bounds_contains(&union_dist, &dist_b));
  }

  #[test]
  fn test_distribution_time_bounds_intersection_contained_by_both() {
    let dist_a = DistributionPlain::range((1.0, 3.0), 1.0);
    let dist_b = DistributionPlain::range((2.0, 4.0), 1.0);

    if let Some((t_min, t_max)) = distribution_time_bounds_intersection(&dist_a, &dist_b) {
      let intersection_dist = DistributionPlain::range((t_min, t_max), 1.0);

      assert!(distribution_time_bounds_contains(&dist_a, &intersection_dist));
      assert!(distribution_time_bounds_contains(&dist_b, &intersection_dist));
    }
  }

  #[test]
  fn test_distribution_time_bounds_negative_ranges() {
    let dist_a = DistributionPlain::range((-5.0, -3.0), 1.0);
    let dist_b = DistributionPlain::range((-4.0, -2.0), 1.0);

    let union_result = distribution_time_bounds_union(&dist_a, &dist_b);
    assert_eq!(Some((-5.0, -2.0)), union_result);

    let intersection_result = distribution_time_bounds_intersection(&dist_a, &dist_b).unwrap();
    assert_eq!((-4.0, -3.0), intersection_result);

    assert!(distribution_time_bounds_overlaps(&dist_a, &dist_b));
  }

  #[test]
  fn test_distribution_time_bounds_zero_width_intervals() {
    let point_a = DistributionPlain::point(2.0, 1.0);
    let point_b = DistributionPlain::point(3.0, 1.0);

    let union_result = distribution_time_bounds_union(&point_a, &point_b);
    assert_eq!(Some((2.0, 3.0)), union_result);

    let intersection_result = distribution_time_bounds_intersection(&point_a, &point_b);
    assert_eq!(None, intersection_result);

    assert!(!distribution_time_bounds_overlaps(&point_a, &point_b));
  }

  #[test]
  fn test_distribution_time_bounds_single_point_overlap() {
    let range_a = DistributionPlain::range((1.0, 2.0), 1.0);
    let range_b = DistributionPlain::range((2.0, 3.0), 1.0);

    let intersection_result = distribution_time_bounds_intersection(&range_a, &range_b).unwrap();
    assert_eq!((2.0, 2.0), intersection_result);

    assert!(distribution_time_bounds_overlaps(&range_a, &range_b));
  }

  #[test]
  fn test_distribution_time_bounds_union_empty_operand_is_identity() {
    let empty = DistributionPlain::empty();
    let range = DistributionPlain::range((1.0, 3.0), 1.0);
    assert_eq!(Some((1.0, 3.0)), distribution_time_bounds_union(&empty, &range));
    assert_eq!(Some((1.0, 3.0)), distribution_time_bounds_union(&range, &empty));
  }

  #[test]
  fn test_distribution_time_bounds_union_both_empty_is_none() {
    let empty = DistributionPlain::empty();
    assert_eq!(None, distribution_time_bounds_union(&empty, &empty));
  }

  #[test]
  fn test_distribution_time_bounds_intersection_empty_operand_is_absorbing() {
    let empty = DistributionPlain::empty();
    let range = DistributionPlain::range((1.0, 3.0), 1.0);
    assert_eq!(None, distribution_time_bounds_intersection(&empty, &range));
    assert_eq!(None, distribution_time_bounds_intersection(&range, &empty));
    assert_eq!(None, distribution_time_bounds_intersection(&empty, &empty));
  }

  #[test]
  fn test_distribution_time_bounds_contains_empty_inner_is_contained() {
    let empty = DistributionPlain::empty();
    let range = DistributionPlain::range((1.0, 3.0), 1.0);
    assert!(distribution_time_bounds_contains(&range, &empty));
    assert!(distribution_time_bounds_contains(&empty, &empty));
  }

  #[test]
  fn test_distribution_time_bounds_contains_empty_outer_contains_no_nonempty() {
    let empty = DistributionPlain::empty();
    let range = DistributionPlain::range((1.0, 3.0), 1.0);
    assert!(!distribution_time_bounds_contains(&empty, &range));
  }

  #[test]
  fn test_distribution_time_bounds_overlaps_empty_never_overlaps() {
    let empty = DistributionPlain::empty();
    let range = DistributionPlain::range((1.0, 3.0), 1.0);
    assert!(!distribution_time_bounds_overlaps(&empty, &range));
    assert!(!distribution_time_bounds_overlaps(&range, &empty));
    assert!(!distribution_time_bounds_overlaps(&empty, &empty));
  }
}
