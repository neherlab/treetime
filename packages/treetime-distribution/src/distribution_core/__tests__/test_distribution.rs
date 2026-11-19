#[cfg(test)]
mod tests {
  use crate::{Distribution, DistributionFunction};
  use eyre::Report;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(   Distribution::Empty,                                                                 None)]
  #[case::point(   Distribution::point(1.0, 5.0),                                                       Some((1.0, 1.0)))]
  #[case::range(   Distribution::range((0.0, 2.0), 3.0),                                                Some((0.0, 2.0)))]
  #[case::function(Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 4.0, 2.0]).unwrap(), Some((0.0, 2.0)))]
  #[trace]
  fn test_distribution_time_bounds_span_the_support(
    #[case] distribution: Distribution,
    #[case] expected: Option<(f64, f64)>,
  ) {
    assert_eq!(expected, distribution.time_bounds());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::peak_and_rounding_above_it(     array![0.0, 1e-15],                    true)]
  #[case::large_values_within_rounding(   array![1000.0, 1000.0 + 1e-13],        true)]
  #[case::zero_probability_at_both_ends(  array![f64::INFINITY, f64::INFINITY],  true)]
  #[case::distinct_probabilities(         array![0.0, 1e-12],                    false)]
  #[case::zero_probability_at_one_end(    array![0.0, f64::INFINITY],            false)]
  #[trace]
  fn test_distribution_function_of_two_equal_probabilities_is_a_range(
    #[case] y: Array1<f64>,
    #[case] is_range: bool,
  ) -> Result<(), Report> {
    let x = array![0.0, 1.0];
    let expected = if is_range {
      Distribution::range((0.0, 1.0), y[1])
    } else {
      Distribution::Function(DistributionFunction::from_range_values((0.0, 1.0), y.clone())?)
    };
    let actual = Distribution::function(x, y)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::at_the_point(               2015.5,                  5.0)]
  #[case::within_rounding_of_a_date(  2015.5 + 1e-12,          5.0)]
  #[case::beyond_rounding_of_a_date(  2015.5 + 1e-11,          f64::INFINITY)]
  #[trace]
  fn test_distribution_point_eval_at_a_date(#[case] t: f64, #[case] expected: f64) -> Result<(), Report> {
    let actual = Distribution::point(2015.5, 5.0).eval(t)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_distribution_point_at_zero_has_zero_probability_at_a_tiny_nonzero_time() -> Result<(), Report> {
    let actual = Distribution::point(0.0, 5.0).eval(1e-17)?;
    assert_eq!(f64::INFINITY, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::hard_sides_end_at_grid(    (BoundaryBehavior::Hard,  BoundaryBehavior::Error),                                        (0.0, 2.0))]
  #[case::soft_sides_are_unbounded(  (BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }), BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 })), (f64::NEG_INFINITY, f64::INFINITY))]
  #[case::hard_approach_ends_at_hard_point(
    (BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: -0.5, b: 1.0 }), BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 3.0, b: 0.0 })),
    (-0.5, 3.0),
  )]
  #[trace]
  fn test_distribution_support_bounds_of_function_end_at_hard_sides(
    #[case] (left, right): (BoundaryBehavior, BoundaryBehavior),
    #[case] expected: (f64, f64),
  ) -> Result<(), Report> {
    let function = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 0.0, 1.0])?
      .with_left_extrap(left)
      .with_right_extrap(right);
    assert_eq!(Some(expected), Distribution::Function(function).support_bounds());
    Ok(())
  }

  #[test]
  fn test_distribution_support_bounds_of_point_range_and_empty_are_time_bounds() {
    assert_eq!(Some((1.0, 1.0)), Distribution::point(1.0, 0.0).support_bounds());
    assert_eq!(Some((0.0, 2.0)), Distribution::range((0.0, 2.0), 0.0).support_bounds());
    assert_eq!(None, Distribution::Empty.support_bounds());
  }

  #[test]
  fn test_distribution_normalize_function_shifts_peak_to_zero() {
    let distribution = Distribution::function(array![0.0, 1.0, 2.0], array![1004.0, 1000.0, 1003.0]).unwrap();
    let Distribution::Function(f) = distribution.normalize().unwrap() else {
      panic!("Expected Function");
    };
    pretty_assert_ulps_eq!(array![4.0, 0.0, 3.0], f.y(), max_ulps = 4);
  }

  #[test]
  fn test_distribution_normalize_is_idempotent() {
    let distribution = Distribution::function(array![0.0, 1.0, 2.0], array![4.0, 0.0, 3.0]).unwrap();
    let once = distribution.normalize().unwrap();
    let twice = once.normalize().unwrap();
    assert_eq!(once, twice);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::point(Distribution::point(2.0, 1000.0),         Distribution::point(2.0, 0.0))]
  #[case::range(Distribution::range((1.0, 3.0), 1000.0),  Distribution::range((1.0, 3.0), 0.0))]
  #[trace]
  fn test_distribution_normalize_constant_maps_to_unit_probability(
    #[case] distribution: Distribution,
    #[case] expected: Distribution,
  ) {
    assert_eq!(expected, distribution.normalize().unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::empty(   Distribution::Empty)]
  #[case::point(   Distribution::point(2.0, f64::INFINITY))]
  #[case::range(   Distribution::range((1.0, 3.0), f64::INFINITY))]
  #[case::function(Distribution::function(array![0.0, 1.0, 2.0], array![f64::INFINITY, f64::INFINITY, f64::INFINITY]).unwrap())]
  #[trace]
  fn test_distribution_normalize_zero_probability_is_empty(#[case] distribution: Distribution) {
    assert_eq!(Distribution::Empty, distribution.normalize().unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::point_nan(    Distribution::point(2.0, f64::NAN),                  "Cannot normalize a distribution point: its peak negative log-likelihood is NaN")]
  #[case::point_neg_inf(Distribution::point(2.0, f64::NEG_INFINITY),         "Cannot normalize a distribution point: its peak negative log-likelihood is -inf")]
  #[case::range_nan(    Distribution::range((1.0, 3.0), f64::NAN),           "Cannot normalize a distribution range: its peak negative log-likelihood is NaN")]
  #[case::function_nan( Distribution::function(array![0.0, 1.0, 2.0], array![1.0, f64::NAN, 3.0]).unwrap(),
                        "Cannot normalize a distribution on [0, 2]: its negative log-likelihood values contain NaN")]
  #[case::function_neg_inf(Distribution::function(array![0.0, 1.0, 2.0], array![1.0, f64::NEG_INFINITY, 3.0]).unwrap(),
                        "Cannot normalize a distribution function on [0, 2]: its peak negative log-likelihood is -inf")]
  #[trace]
  fn test_distribution_normalize_rejects_invalid_likelihood(#[case] distribution: Distribution, #[case] expected: &str) {
    assert_error!(distribution.normalize(), expected);
  }
}
