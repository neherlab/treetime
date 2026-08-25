#[cfg(test)]
mod tests {
  use crate::coalescent::population_size::effective_population_size;
  use crate::pretty_assert_ulps_eq;
  use rstest::rstest;

  // Oracle: the definition N_e = Tc * gen_per_year. Expected values are the products worked out by
  // hand, not recomputed with the same expression the function uses.
  #[rustfmt::skip]
  #[rstest]
  #[case::default_factor(   2.0, 50.0, 100.0)]
  #[case::skyline_segment(  3.0, 50.0, 150.0)]
  #[case::small_timescale(  0.5,  2.0,   1.0)]
  #[case::fractional_factor(4.0,  0.25,  1.0)]
  #[trace]
  fn test_population_size_value_equals_tc_times_factor(
    #[case] tc: f64,
    #[case] gen_per_year: f64,
    #[case] expected: f64,
  ) {
    pretty_assert_ulps_eq!(expected, effective_population_size(tc, gen_per_year), max_ulps = 4);
  }

  // A confidence band on Tc maps to a band on N_e by scaling each bound with the same factor, so
  // the band still brackets the point estimate and its width scales linearly. Expected values are
  // the hand-computed products for gen_per_year = 50.0.
  #[test]
  fn test_population_size_band_scales_linearly() {
    let gen_per_year = 50.0;
    let (tc_lower, tc_value, tc_upper) = (1.5, 3.0, 6.0);

    let ne_lower = effective_population_size(tc_lower, gen_per_year);
    let ne_value = effective_population_size(tc_value, gen_per_year);
    let ne_upper = effective_population_size(tc_upper, gen_per_year);

    pretty_assert_ulps_eq!(75.0, ne_lower, max_ulps = 4);
    pretty_assert_ulps_eq!(150.0, ne_value, max_ulps = 4);
    pretty_assert_ulps_eq!(300.0, ne_upper, max_ulps = 4);

    // Band brackets the point estimate and its width is the Tc width scaled by the factor.
    assert!(ne_lower <= ne_value && ne_value <= ne_upper);
    pretty_assert_ulps_eq!(gen_per_year * (tc_upper - tc_lower), ne_upper - ne_lower, max_ulps = 4);
  }

  // Fixed Tc has no confidence band, so the report carries only the point conversion. Here the
  // default gen_per_year = 50.0 rescales a user-supplied Tc into a single N_e value.
  #[test]
  fn test_population_size_fixed_mode_value_only() {
    pretty_assert_ulps_eq!(60.0, effective_population_size(1.2, 50.0), max_ulps = 4);
  }
}
