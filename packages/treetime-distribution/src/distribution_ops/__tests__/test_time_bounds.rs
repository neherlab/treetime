#[cfg(test)]
mod tests {
  use crate::distribution_ops::time_bounds::distribution_support_n_points;
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
}
