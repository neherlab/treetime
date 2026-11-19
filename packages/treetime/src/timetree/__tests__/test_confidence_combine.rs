#[cfg(test)]
mod tests {
  use crate::timetree::confidence::{combine_in_quadrature, confidence_bounds};
  use crate::timetree::max_posterior_region::MaxPosteriorRegion;
  use rstest::rstest;
  use treetime_utils::pretty_assert_relative_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::no_contribution_keeps_the_limits(  (5.0, 15.0), None,              None,              (5.0, 15.0))]
  #[case::first_contribution_alone(          (0.0, 20.0), Some((8.0, 12.0)), None,              (8.0, 12.0))]
  #[case::second_contribution_alone(         (0.0, 20.0), None,              Some((7.0, 13.0)), (7.0, 13.0))]
  #[case::single_contribution_clipped(       (9.0, 11.0), Some((8.0, 12.0)), None,              (9.0, 11.0))]
  #[case::combined_contributions_clipped(    (5.0, 15.0), Some((0.0, 20.0)), Some((0.0, 20.0)), (5.0, 15.0))]
  #[trace]
  fn test_combine_in_quadrature_takes_a_single_contribution_and_clips_to_the_limits(
    #[case] limits: (f64, f64),
    #[case] c1: Option<(f64, f64)>,
    #[case] c2: Option<(f64, f64)>,
    #[case] expected: (f64, f64),
  ) {
    let actual = combine_in_quadrature(10.0, limits, c1, c2);
    pretty_assert_relative_eq!(expected.0, actual.0);
    pretty_assert_relative_eq!(expected.1, actual.1);
  }

  #[test]
  fn test_combine_in_quadrature_adds_the_deviations_of_both_contributions_in_quadrature() {
    let actual = combine_in_quadrature(10.0, (0.0, 20.0), Some((8.0, 12.0)), Some((7.0, 13.0)));
    let expected_dev = (4.0_f64 + 9.0).sqrt();
    pretty_assert_relative_eq!(10.0 - expected_dev, actual.0, epsilon = 1e-10);
    pretty_assert_relative_eq!(10.0 + expected_dev, actual.1, epsilon = 1e-10);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::without_rate_uncertainty(None)]
  #[case::with_rate_uncertainty(   Some([2019.0, 2020.5, 2022.0]))]
  #[case::with_shifted_rate_dates( Some([2020.0, 2020.2, 2020.4]))]
  #[trace]
  fn test_confidence_bounds_point_region_is_the_date_without_limits(#[case] rate_dates: Option<[f64; 3]>) {
    let date = 2020.5;
    let unbounded = (f64::NEG_INFINITY, f64::INFINITY);

    let interval = confidence_bounds(date, unbounded, Some(MaxPosteriorRegion::Point), rate_dates.as_ref());

    let v0_delta_region = (date, date);
    assert_eq!(v0_delta_region, interval);
  }
}
