#[cfg(test)]
mod tests {
  use crate::distribution_ops::time_bounds::distribution_support_n_points;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_grid::MaxGridPoints;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::fractional_ceils_up(     (0.0, 2.4),   1.0, 4)]
  #[case::small_fraction_ceils_up( (0.0, 2.1),   1.0, 4)]
  #[case::exact_multiple_no_extra( (0.0, 3.0),   1.0, 4)]
  #[case::minimum_two(             (0.0, 0.4),   1.0, 2)]
  #[case::exactly_at_the_limit(    (0.0, 999.0), 1.0, 1_000)]
  #[trace]
  fn test_distribution_support_n_points_uses_spacing_contract(
    #[case] bounds: (f64, f64),
    #[case] dx: f64,
    #[case] expected: usize,
  ) {
    let actual = distribution_support_n_points(bounds, dx, MaxGridPoints::new(1_000).unwrap()).unwrap();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_point_above_the_limit((0.0, 1000.0), 1.0,    "A grid over [0, 1000] with spacing 1 needs 1001 points, more than the limit of 1000")]
  #[case::ratio_overflows_usize(    (0.0, 1.0),    1e-320, "A grid over [0, 1] with spacing 1.0e-320 needs more points than the limit of 1000")]
  #[trace]
  fn test_distribution_support_n_points_stops_above_the_limit(
    #[case] bounds: (f64, f64),
    #[case] dx: f64,
    #[case] expected: &str,
  ) {
    assert_error!(
      distribution_support_n_points(bounds, dx, MaxGridPoints::new(1_000).unwrap()),
      expected
    );
  }
}
