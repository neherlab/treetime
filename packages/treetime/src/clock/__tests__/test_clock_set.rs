#[cfg(test)]
mod tests {
  use crate::clock::clock_set::ClockSet;

  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rstest]
  #[case::zero_branch(0.0)]
  #[case::short_branch(0.1)]
  #[case::long_branch(10.0)]
  #[trace]
  fn test_clock_set_dateless_leaf_contribution_is_zero(#[case] branch_length: f64) {
    let expected = ClockSet::outlier_contribution();
    let actual = ClockSet::leaf_contribution(None).shift(branch_length);

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_clock_set_shift_equals_sums_of_lengthened_divergences() {
    let tips = [(2019.0, 0.01), (2020.5, 0.03), (2021.0, 0.05)];
    let l = 0.02;
    let set = |shift: f64| {
      tips.iter().fold(ClockSet::default(), |acc, &(t, d)| {
        acc + ClockSet::leaf_contribution(Some(t)).shift(d + shift)
      })
    };

    let shifted = set(0.0).shift(l);
    let direct = set(l);
    pretty_assert_ulps_eq!(shifted.d_sum(), direct.d_sum(), max_ulps = 4);
    pretty_assert_ulps_eq!(shifted.dt_sum(), direct.dt_sum(), max_ulps = 4);
    pretty_assert_ulps_eq!(shifted.dsq_sum(), direct.dsq_sum(), max_ulps = 4);
    assert_eq!(
      (shifted.t_sum(), shifted.tsq_sum(), shifted.norm()),
      (direct.t_sum(), direct.tsq_sum(), direct.norm())
    );
  }

  #[test]
  fn test_clock_set_chisq_fixed_rate_matches_zero_rate_formula() {
    let cs = ClockSet::leaf_contribution(Some(0.0)).shift(1.0) + ClockSet::leaf_contribution(Some(1.0)).shift(3.0);

    // residuals -1, +1 about the mean divergence 2: RSS / 2 = 1
    pretty_assert_ulps_eq!(cs.chisq_fixed_rate(0.0), 1.0, max_ulps = 4);
  }

  #[test]
  fn test_clock_set_chisq_fixed_rate_matches_nonzero_rate_formula() {
    let cs = ClockSet::leaf_contribution(Some(0.0)).shift(1.0) + ClockSet::leaf_contribution(Some(1.0)).shift(3.0);

    // d - t = 1, 2 about their mean 1.5: RSS / 2 = 0.25
    pretty_assert_ulps_eq!(cs.chisq_fixed_rate(1.0), 0.25, max_ulps = 4);
  }
}
