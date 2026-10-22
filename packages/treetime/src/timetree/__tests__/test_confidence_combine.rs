#[cfg(test)]
mod tests {
  use crate::timetree::confidence::{combine_in_quadrature, confidence_bounds};
  use crate::timetree::max_posterior_region::MaxPosteriorRegion;
  use rstest::rstest;
  use treetime_utils::pretty_assert_relative_eq;

  #[test]
  fn test_combine_in_quadrature_no_contributions() {
    let result = combine_in_quadrature(10.0, (5.0, 15.0), None, None);
    pretty_assert_relative_eq!(result.0, 5.0);
    pretty_assert_relative_eq!(result.1, 15.0);
  }

  #[test]
  fn test_combine_in_quadrature_single_contribution() {
    let result = combine_in_quadrature(10.0, (0.0, 20.0), Some((8.0, 12.0)), None);
    pretty_assert_relative_eq!(result.0, 8.0);
    pretty_assert_relative_eq!(result.1, 12.0);

    let result = combine_in_quadrature(10.0, (0.0, 20.0), None, Some((7.0, 13.0)));
    pretty_assert_relative_eq!(result.0, 7.0);
    pretty_assert_relative_eq!(result.1, 13.0);
  }

  #[test]
  fn test_combine_in_quadrature_quadrature() {
    let result = combine_in_quadrature(10.0, (0.0, 20.0), Some((8.0, 12.0)), Some((7.0, 13.0)));
    let expected_dev = (4.0_f64 + 9.0).sqrt();
    pretty_assert_relative_eq!(result.0, 10.0 - expected_dev, epsilon = 1e-10);
    pretty_assert_relative_eq!(result.1, 10.0 + expected_dev, epsilon = 1e-10);
  }

  #[test]
  fn test_combine_in_quadrature_clipped_to_limits() {
    let result = combine_in_quadrature(10.0, (5.0, 15.0), Some((0.0, 20.0)), Some((0.0, 20.0)));
    pretty_assert_relative_eq!(result.0, 5.0);
    pretty_assert_relative_eq!(result.1, 15.0);
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
