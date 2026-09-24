#[cfg(test)]
mod tests {
  use crate::Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::policy::{NegLog, Plain};
  use ndarray::array;
  use rstest::rstest;
  use treetime_utils::{assert_error, make_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  #[test]
  fn test_distribution_time_bounds_empty_is_none() {
    assert_eq!(None, Distribution::<Plain>::Empty.time_bounds());
  }

  #[test]
  fn test_distribution_time_bounds_point_is_degenerate_interval() {
    assert_eq!(Some((1.0, 1.0)), Distribution::<Plain>::point(1.0, 5.0).time_bounds());
  }

  #[test]
  fn test_distribution_time_bounds_range_spans_endpoints() {
    assert_eq!(
      Some((0.0, 2.0)),
      Distribution::<Plain>::range((0.0, 2.0), 3.0).time_bounds()
    );
  }

  #[test]
  fn test_distribution_time_bounds_function_spans_grid() {
    let func = Distribution::<Plain>::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]).unwrap();
    assert_eq!(Some((0.0, 2.0)), func.time_bounds());
  }

  #[test]
  fn test_distribution_y_formula_evaluates_bounds() {
    let formula = Distribution::<NegLog>::Formula(DistributionFormula::new(|t| Ok(2.0 * t), 1.0, 3.0));
    assert_eq!(array![2.0, 6.0], formula.y().unwrap());
  }

  #[test]
  fn test_distribution_y_formula_propagates_evaluation_error() {
    let formula =
      Distribution::<NegLog>::Formula(DistributionFormula::new(|t| make_error!("no value at {t}"), 1.0, 3.0));
    assert_error!(
      formula.y(),
      "When evaluating a formula distribution at its bounds [1, 3]: no value at 1"
    );
  }

  #[test]
  fn test_distribution_neglog_normalize_function_shifts_peak_to_zero() {
    let distribution = Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![1004.0, 1000.0, 1003.0]).unwrap();
    let Distribution::Function(f) = distribution.normalize().unwrap() else {
      panic!("Expected Function");
    };
    pretty_assert_ulps_eq!(array![4.0, 0.0, 3.0], f.y(), max_ulps = 4);
  }

  #[test]
  fn test_distribution_neglog_normalize_is_idempotent() {
    let distribution = Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![4.0, 0.0, 3.0]).unwrap();
    let once = distribution.normalize().unwrap();
    let twice = once.normalize().unwrap();
    assert_eq!(once, twice);
  }

  #[test]
  fn test_distribution_neglog_normalize_point_maps_to_unit_probability() {
    let expected = Distribution::<NegLog>::point(2.0, 0.0);
    assert_eq!(
      expected,
      Distribution::<NegLog>::point(2.0, 1000.0).normalize().unwrap()
    );
  }

  #[test]
  fn test_distribution_neglog_normalize_range_maps_to_unit_probability() {
    let expected = Distribution::<NegLog>::range((1.0, 3.0), 0.0);
    assert_eq!(
      expected,
      Distribution::<NegLog>::range((1.0, 3.0), 1000.0).normalize().unwrap()
    );
  }

  #[test]
  fn test_distribution_neglog_normalize_formula_discretizes_and_shifts_peak() {
    let formula = Distribution::<NegLog>::Formula(DistributionFormula::new(|t| Ok((t - 1.0).powi(2) + 7.0), 0.0, 2.0));
    let Distribution::Function(f) = formula.normalize().unwrap() else {
      panic!("Expected Function");
    };
    let grid_min_offset = (1.0_f64 / 199.0).powi(2);
    let expected = f.t().mapv(|t| (t - 1.0).powi(2) - grid_min_offset);
    pretty_assert_abs_diff_eq!(expected, f.y(), epsilon = 1e-12);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(   Distribution::<NegLog>::Empty)]
  #[case::point(   Distribution::<NegLog>::point(2.0, f64::INFINITY))]
  #[case::range(   Distribution::<NegLog>::range((1.0, 3.0), f64::INFINITY))]
  #[case::function(Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![f64::INFINITY, f64::INFINITY, f64::INFINITY]).unwrap())]
  #[trace]
  fn test_distribution_neglog_normalize_zero_probability_is_empty(#[case] distribution: Distribution<NegLog>) {
    assert_eq!(Distribution::<NegLog>::Empty, distribution.normalize().unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::point_nan(    Distribution::<NegLog>::point(2.0, f64::NAN),                  "Cannot normalize a distribution point: its peak negative log-likelihood is NaN")]
  #[case::point_neg_inf(Distribution::<NegLog>::point(2.0, f64::NEG_INFINITY),         "Cannot normalize a distribution point: its peak negative log-likelihood is -inf")]
  #[case::range_nan(    Distribution::<NegLog>::range((1.0, 3.0), f64::NAN),           "Cannot normalize a distribution range: its peak negative log-likelihood is NaN")]
  #[case::function_nan( Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![1.0, f64::NAN, 3.0]).unwrap(),
                        "Cannot normalize a distribution on [0, 2]: its negative log-likelihood values contain NaN")]
  #[case::function_neg_inf(Distribution::<NegLog>::function(array![0.0, 1.0, 2.0], array![1.0, f64::NEG_INFINITY, 3.0]).unwrap(),
                        "Cannot normalize a distribution function on [0, 2]: its peak negative log-likelihood is -inf")]
  #[trace]
  fn test_distribution_neglog_normalize_rejects_invalid_likelihood(#[case] distribution: Distribution<NegLog>, #[case] expected: &str) {
    assert_error!(distribution.normalize(), expected);
  }

  #[test]
  fn test_distribution_neglog_normalize_formula_propagates_evaluation_error() {
    let formula =
      Distribution::<NegLog>::Formula(DistributionFormula::new(|t| make_error!("no value at {t}"), 1.0, 3.0));
    assert_error!(
      formula.normalize(),
      "When discretizing a formula distribution on [1, 3] for normalization: no value at 1"
    );
  }

  #[test]
  fn test_distribution_likely_time_formula_selects_min_neglog() {
    let formula = Distribution::<NegLog>::Formula(DistributionFormula::new(|t| Ok((t - 1.0).abs()), 0.0, 1.99));
    pretty_assert_ulps_eq!(1.0, formula.likely_time().unwrap().unwrap(), max_ulps = 4);
  }

  #[test]
  fn test_distribution_likely_time_formula_propagates_evaluation_error() {
    let formula =
      Distribution::<NegLog>::Formula(DistributionFormula::new(|t| make_error!("no value at {t}"), 1.0, 3.0));
    assert_error!(
      formula.likely_time(),
      "When finding the most likely time of a formula distribution on [1, 3]: no value at 1"
    );
  }

  #[test]
  fn test_distribution_likely_time_formula_rejects_nan() {
    let formula = Distribution::<NegLog>::Formula(DistributionFormula::new(|_| Ok(f64::NAN), 1.0, 3.0));
    assert_error!(
      formula.likely_time(),
      "Cannot find the most likely time of a formula distribution on [1, 3]: its values contain NaN"
    );
  }
}
