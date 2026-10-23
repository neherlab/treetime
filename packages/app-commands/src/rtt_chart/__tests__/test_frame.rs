#[cfg(test)]
mod tests {
  use crate::rtt_chart::__tests__::test_support::tests::{estimated_model, fitted_tips, tip};
  use crate::rtt_chart::data::RttChartData;
  use crate::rtt_chart::frame::{RootPlacement, RttFrame};
  use pretty_assertions::assert_eq;

  #[test]
  fn test_rtt_frame_extrapolates_to_a_root_close_to_the_samples() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.5), (2003.0, 1.5), (2005.0, 2.5)])).unwrap();
    let expected = RttFrame {
      dates: (2000.0 - 5.0 * 0.04, 2005.0 + 5.0 * 0.04),
      divergences: (0.0 - 2.5 * 0.06, 2.5 + 2.5 * 0.06),
      fit: [(2001.0, 0.5), (2005.0, 2.5)],
      extrapolation: Some([(2000.0, 0.0), (2001.0, 0.5)]),
      root: RootPlacement::OnChart(2000.0),
    };
    assert_eq!(expected, RttFrame::new(&data, &estimated_model(0.5, -1000.0, 0.9)));
  }

  #[test]
  fn test_rtt_frame_points_toward_a_root_far_before_the_samples() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.5), (2003.0, 1.5)])).unwrap();
    let frame = RttFrame::new(&data, &estimated_model(0.5, -1000.5 + 0.5 * 4.0, 0.9));
    let left = 2001.0 - 2.0 * 0.04;
    let expected = (
      (left, 2003.0 + 2.0 * 0.04),
      Some([(left, 0.5 * left - 998.5), (2001.0, 2001.0 * 0.5 - 998.5)]),
      RootPlacement::BeforeChart(998.5 / 0.5),
    );
    assert_eq!(expected, (frame.dates, frame.extrapolation, frame.root));
  }

  #[test]
  fn test_rtt_frame_without_positive_rate_has_no_root() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 2.5), (2005.0, 0.5)])).unwrap();
    let frame = RttFrame::new(&data, &estimated_model(-0.5, 1003.0, -0.9));
    assert_eq!((None, RootPlacement::Undefined), (frame.extrapolation, frame.root));
  }

  #[test]
  fn test_rtt_frame_with_the_root_inside_the_fitted_range_draws_no_extrapolation() {
    let data = RttChartData::new(&fitted_tips(&[(1999.0, 0.0), (2001.0, 0.5), (2005.0, 2.5)])).unwrap();
    let frame = RttFrame::new(&data, &estimated_model(0.5, -1000.0, 0.9));
    assert_eq!(
      (None, RootPlacement::OnChart(2000.0)),
      (frame.extrapolation, frame.root)
    );
  }

  #[test]
  fn test_rtt_frame_covers_excluded_samples() {
    let data = RttChartData::new(&[
      tip(Some(2001.0), 0.5, false),
      tip(Some(2005.0), 2.5, false),
      tip(Some(2006.0), 9.5, true),
    ])
    .unwrap();
    let frame = RttFrame::new(&data, &estimated_model(0.5, -1000.0, 0.9));
    assert_eq!(
      (
        (2000.0 - 6.0 * 0.04, 2006.0 + 6.0 * 0.04),
        (0.0 - 9.5 * 0.06, 9.5 + 9.5 * 0.06)
      ),
      (frame.dates, frame.divergences)
    );
  }
}
