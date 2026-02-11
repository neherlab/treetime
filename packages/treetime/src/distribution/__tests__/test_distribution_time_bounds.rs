#[cfg(test)]
mod tests {
  use crate::distribution::distribution::DistributionPlain as Distribution;
  use crate::distribution::distribution_time_bounds::{
    distribution_time_bounds_contains, distribution_time_bounds_intersection, distribution_time_bounds_overlaps,
    distribution_time_bounds_union,
  };
  use ndarray::array;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::two_points_ordered(Distribution::point(1.0, 1.0),           Distribution::point(2.0, 1.0),          (1.0, 2.0)  )]
  #[case::two_points_reversed(Distribution::point(2.0, 1.0),           Distribution::point(1.0, 1.0),          (1.0, 2.0)  )]
  #[case::overlapping_ranges(Distribution::range((1.0, 3.0), 1.0),    Distribution::range((2.0, 4.0), 1.0),   (1.0, 4.0)  )]
  #[case::overlapping_ranges_reversed(Distribution::range((2.0, 4.0), 1.0),    Distribution::range((1.0, 3.0), 1.0),   (1.0, 4.0)  )]
  #[case::disjoint_ranges(Distribution::range((1.0, 2.0), 1.0),    Distribution::range((3.0, 4.0), 1.0),   (1.0, 4.0)  )]
  #[case::nested_ranges(Distribution::range((1.0, 5.0), 1.0),    Distribution::range((2.0, 3.0), 1.0),   (1.0, 5.0)  )]
  #[case::point_and_range(Distribution::point(1.0, 1.0),           Distribution::range((2.0, 3.0), 1.0),   (1.0, 3.0)  )]
  #[case::two_functions(
    Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
    Distribution::function(array![1.5, 2.5, 3.5], array![1.0, 2.0, 3.0]).unwrap(),
    (0.0, 3.5)
  )]
  #[case::same_point(Distribution::point(5.0, 1.0),           Distribution::point(5.0, 1.0),          (5.0, 5.0)  )]
  #[case::negative_and_positive(Distribution::range((-2.0, -1.0), 1.0),  Distribution::range((1.0, 2.0), 1.0),   (-2.0, 2.0) )]
  #[trace]
  fn test_distribution_time_bounds_union(
    #[case] dist_a: Distribution,
    #[case] dist_b: Distribution,
    #[case] expected: (f64, f64),
  ) {
    let actual = distribution_time_bounds_union(&dist_a, &dist_b);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::overlapping_ranges(Distribution::range((1.0, 3.0), 1.0),   Distribution::range((2.0, 4.0), 1.0),   Some((2.0, 3.0)))]
  #[case::overlapping_ranges_reversed(Distribution::range((2.0, 4.0), 1.0),   Distribution::range((1.0, 3.0), 1.0),   Some((2.0, 3.0)))]
  #[case::disjoint_ranges(Distribution::range((1.0, 2.0), 1.0),   Distribution::range((3.0, 4.0), 1.0),   None            )]
  #[case::disjoint_ranges_reversed(Distribution::range((3.0, 4.0), 1.0),   Distribution::range((1.0, 2.0), 1.0),   None            )]
  #[case::nested_outer_first(Distribution::range((1.0, 5.0), 1.0),   Distribution::range((2.0, 3.0), 1.0),   Some((2.0, 3.0)))]
  #[case::nested_inner_first(Distribution::range((2.0, 3.0), 1.0),   Distribution::range((1.0, 5.0), 1.0),   Some((2.0, 3.0)))]
  #[case::point_inside_range(Distribution::point(2.5, 1.0),          Distribution::range((2.0, 3.0), 1.0),   Some((2.5, 2.5)))]
  #[case::range_contains_point(Distribution::range((2.0, 3.0), 1.0),   Distribution::point(2.5, 1.0),          Some((2.5, 2.5)))]
  #[case::point_outside_range(Distribution::point(1.0, 1.0),          Distribution::range((2.0, 3.0), 1.0),   None            )]
  #[case::two_functions(
    Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap(),
    Distribution::function(array![1.5, 2.5, 3.5], array![1.0, 2.0, 3.0]).unwrap(),
    Some((1.5, 2.0))
  )]
  #[case::adjacent_ranges(Distribution::range((1.0, 2.0), 1.0),   Distribution::range((2.0, 3.0), 1.0),   Some((2.0, 2.0)))]
  #[case::same_point(Distribution::point(5.0, 1.0),          Distribution::point(5.0, 1.0),          Some((5.0, 5.0)))]
  #[case::different_points(Distribution::point(5.0, 1.0),          Distribution::point(6.0, 1.0),          None            )]
  #[trace]
  fn test_distribution_time_bounds_intersection(
    #[case] dist_a: Distribution,
    #[case] dist_b: Distribution,
    #[case] expected: Option<(f64, f64)>,
  ) {
    let actual = distribution_time_bounds_intersection(&dist_a, &dist_b);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::outer_contains_inner(Distribution::range((1.0, 5.0), 1.0),   Distribution::range((2.0, 3.0), 1.0),   true )]
  #[case::inner_not_contains_outer(Distribution::range((2.0, 3.0), 1.0),   Distribution::range((1.0, 5.0), 1.0),   false)]
  #[case::range_contains_point(Distribution::range((1.0, 5.0), 1.0),   Distribution::point(3.0, 1.0),          true )]
  #[case::point_not_contains_range(Distribution::point(3.0, 1.0),          Distribution::range((1.0, 5.0), 1.0),   false)]
  #[case::overlapping_not_contained(Distribution::range((1.0, 3.0), 1.0),   Distribution::range((2.0, 4.0), 1.0),   false)]
  #[case::same_range(Distribution::range((1.0, 5.0), 1.0),   Distribution::range((1.0, 5.0), 1.0),   true )]
  #[case::same_point(Distribution::point(5.0, 1.0),          Distribution::point(5.0, 1.0),          true )]
  #[case::inner_extends_left(Distribution::range((1.0, 5.0), 1.0),   Distribution::range((0.5, 3.0), 1.0),   false)]
  #[case::inner_extends_right(Distribution::range((1.0, 5.0), 1.0),   Distribution::range((3.0, 5.5), 1.0),   false)]
  #[case::function_contains_function(
    Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0], array![1.0, 2.0, 3.0, 2.0, 1.0, 0.5]).unwrap(),
    Distribution::function(array![1.0, 2.0, 3.0], array![1.0, 2.0, 3.0]).unwrap(),
    true
  )]
  #[trace]
  fn test_distribution_time_bounds_contains(
    #[case] outer: Distribution,
    #[case] inner: Distribution,
    #[case] expected: bool,
  ) {
    let actual = distribution_time_bounds_contains(&outer, &inner);
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::overlapping_ranges(Distribution::range((1.0, 3.0), 1.0), Distribution::range((2.0, 4.0), 1.0), true )]
  #[case::disjoint_ranges(Distribution::range((1.0, 2.0), 1.0), Distribution::range((3.0, 4.0), 1.0), false)]
  #[case::adjacent_ranges(Distribution::range((1.0, 2.0), 1.0), Distribution::range((2.0, 3.0), 1.0), true )]
  #[case::nested_ranges(Distribution::range((1.0, 5.0), 1.0), Distribution::range((2.0, 3.0), 1.0), true )]
  #[case::point_inside_range(Distribution::point(2.5, 1.0),        Distribution::range((2.0, 3.0), 1.0), true )]
  #[case::point_outside_range(Distribution::point(1.0, 1.0),        Distribution::range((2.0, 3.0), 1.0), false)]
  #[case::same_point(Distribution::point(5.0, 1.0),        Distribution::point(5.0, 1.0),        true )]
  #[case::different_points(Distribution::point(5.0, 1.0),        Distribution::point(6.0, 1.0),        false)]
  #[trace]
  fn test_distribution_time_bounds_overlaps(
    #[case] dist_a: Distribution,
    #[case] dist_b: Distribution,
    #[case] expected: bool,
  ) {
    let actual = distribution_time_bounds_overlaps(&dist_a, &dist_b);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_distribution_time_bounds_union_commutativity() {
    let dist_a = Distribution::range((1.0, 3.0), 1.0);
    let dist_b = Distribution::range((2.0, 4.0), 1.0);

    let result_ab = distribution_time_bounds_union(&dist_a, &dist_b);
    let result_ba = distribution_time_bounds_union(&dist_b, &dist_a);

    assert_eq!(result_ab, result_ba);
  }

  #[test]
  fn test_distribution_time_bounds_intersection_commutativity() {
    let dist_a = Distribution::range((1.0, 3.0), 1.0);
    let dist_b = Distribution::range((2.0, 4.0), 1.0);

    let result_ab = distribution_time_bounds_intersection(&dist_a, &dist_b);
    let result_ba = distribution_time_bounds_intersection(&dist_b, &dist_a);

    assert_eq!(result_ab, result_ba);
  }

  #[test]
  fn test_distribution_time_bounds_overlaps_symmetry() {
    let dist_a = Distribution::range((1.0, 3.0), 1.0);
    let dist_b = Distribution::range((2.0, 4.0), 1.0);

    let result_ab = distribution_time_bounds_overlaps(&dist_a, &dist_b);
    let result_ba = distribution_time_bounds_overlaps(&dist_b, &dist_a);

    assert_eq!(result_ab, result_ba);
  }

  #[test]
  fn test_distribution_time_bounds_union_associativity() {
    let dist_a = Distribution::range((1.0, 2.0), 1.0);
    let dist_b = Distribution::range((3.0, 4.0), 1.0);
    let dist_c = Distribution::range((5.0, 6.0), 1.0);

    let (t_min_ab, t_max_ab) = distribution_time_bounds_union(&dist_a, &dist_b);
    let dist_ab = Distribution::range((t_min_ab, t_max_ab), 1.0);
    let result_ab_c = distribution_time_bounds_union(&dist_ab, &dist_c);

    let (t_min_bc, t_max_bc) = distribution_time_bounds_union(&dist_b, &dist_c);
    let dist_bc = Distribution::range((t_min_bc, t_max_bc), 1.0);
    let result_a_bc = distribution_time_bounds_union(&dist_a, &dist_bc);

    assert_eq!(result_ab_c, result_a_bc);
  }

  #[test]
  fn test_distribution_time_bounds_intersection_associativity() {
    let dist_a = Distribution::range((1.0, 6.0), 1.0);
    let dist_b = Distribution::range((2.0, 5.0), 1.0);
    let dist_c = Distribution::range((3.0, 4.0), 1.0);

    let int_ab = distribution_time_bounds_intersection(&dist_a, &dist_b).unwrap();
    let dist_ab = Distribution::range(int_ab, 1.0);
    let result_ab_c = distribution_time_bounds_intersection(&dist_ab, &dist_c);

    let int_bc = distribution_time_bounds_intersection(&dist_b, &dist_c).unwrap();
    let dist_bc = Distribution::range(int_bc, 1.0);
    let result_a_bc = distribution_time_bounds_intersection(&dist_a, &dist_bc);

    assert_eq!(result_ab_c, result_a_bc);
  }

  #[test]
  fn test_distribution_time_bounds_idempotence() {
    let dist = Distribution::range((1.0, 3.0), 1.0);

    let union_result = distribution_time_bounds_union(&dist, &dist);
    let expected = dist.time_bounds();
    assert_eq!(expected, union_result);

    let intersection_result = distribution_time_bounds_intersection(&dist, &dist).unwrap();
    assert_eq!(expected, intersection_result);
  }

  #[test]
  fn test_distribution_time_bounds_contains_reflexivity() {
    let dist = Distribution::range((1.0, 3.0), 1.0);
    let actual = distribution_time_bounds_contains(&dist, &dist);
    assert!(actual);
  }

  #[test]
  fn test_distribution_time_bounds_contains_transitivity() {
    let dist_outer = Distribution::range((1.0, 10.0), 1.0);
    let dist_middle = Distribution::range((2.0, 8.0), 1.0);
    let dist_inner = Distribution::range((3.0, 7.0), 1.0);

    assert!(distribution_time_bounds_contains(&dist_outer, &dist_middle));
    assert!(distribution_time_bounds_contains(&dist_middle, &dist_inner));
    assert!(distribution_time_bounds_contains(&dist_outer, &dist_inner));
  }

  #[test]
  fn test_distribution_time_bounds_overlaps_reflexivity() {
    let dist = Distribution::range((1.0, 3.0), 1.0);
    let actual = distribution_time_bounds_overlaps(&dist, &dist);
    assert!(actual);
  }

  #[test]
  fn test_distribution_time_bounds_intersection_implies_overlap() {
    let dist_a = Distribution::range((1.0, 3.0), 1.0);
    let dist_b = Distribution::range((2.0, 4.0), 1.0);

    let has_intersection = distribution_time_bounds_intersection(&dist_a, &dist_b).is_some();
    let overlaps = distribution_time_bounds_overlaps(&dist_a, &dist_b);

    assert_eq!(has_intersection, overlaps);
  }

  #[test]
  fn test_distribution_time_bounds_contains_implies_overlap() {
    let outer = Distribution::range((1.0, 5.0), 1.0);
    let inner = Distribution::range((2.0, 3.0), 1.0);

    assert!(distribution_time_bounds_contains(&outer, &inner));
    assert!(distribution_time_bounds_overlaps(&outer, &inner));
  }

  #[test]
  fn test_distribution_time_bounds_union_contains_both() {
    let dist_a = Distribution::range((1.0, 3.0), 1.0);
    let dist_b = Distribution::range((2.0, 4.0), 1.0);

    let (t_min, t_max) = distribution_time_bounds_union(&dist_a, &dist_b);
    let union_dist = Distribution::range((t_min, t_max), 1.0);

    assert!(distribution_time_bounds_contains(&union_dist, &dist_a));
    assert!(distribution_time_bounds_contains(&union_dist, &dist_b));
  }

  #[test]
  fn test_distribution_time_bounds_intersection_contained_by_both() {
    let dist_a = Distribution::range((1.0, 3.0), 1.0);
    let dist_b = Distribution::range((2.0, 4.0), 1.0);

    if let Some((t_min, t_max)) = distribution_time_bounds_intersection(&dist_a, &dist_b) {
      let intersection_dist = Distribution::range((t_min, t_max), 1.0);

      assert!(distribution_time_bounds_contains(&dist_a, &intersection_dist));
      assert!(distribution_time_bounds_contains(&dist_b, &intersection_dist));
    }
  }

  #[test]
  fn test_distribution_time_bounds_negative_ranges() {
    let dist_a = Distribution::range((-5.0, -3.0), 1.0);
    let dist_b = Distribution::range((-4.0, -2.0), 1.0);

    let union_result = distribution_time_bounds_union(&dist_a, &dist_b);
    assert_eq!((-5.0, -2.0), union_result);

    let intersection_result = distribution_time_bounds_intersection(&dist_a, &dist_b).unwrap();
    assert_eq!((-4.0, -3.0), intersection_result);

    assert!(distribution_time_bounds_overlaps(&dist_a, &dist_b));
  }

  #[test]
  fn test_distribution_time_bounds_zero_width_intervals() {
    let point_a = Distribution::point(2.0, 1.0);
    let point_b = Distribution::point(3.0, 1.0);

    let union_result = distribution_time_bounds_union(&point_a, &point_b);
    assert_eq!((2.0, 3.0), union_result);

    let intersection_result = distribution_time_bounds_intersection(&point_a, &point_b);
    assert_eq!(None, intersection_result);

    assert!(!distribution_time_bounds_overlaps(&point_a, &point_b));
  }

  #[test]
  fn test_distribution_time_bounds_single_point_overlap() {
    let range_a = Distribution::range((1.0, 2.0), 1.0);
    let range_b = Distribution::range((2.0, 3.0), 1.0);

    let intersection_result = distribution_time_bounds_intersection(&range_a, &range_b).unwrap();
    assert_eq!((2.0, 2.0), intersection_result);

    assert!(distribution_time_bounds_overlaps(&range_a, &range_b));
  }
}
