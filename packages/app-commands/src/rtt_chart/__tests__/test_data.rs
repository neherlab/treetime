#[cfg(test)]
mod tests {
  use crate::rtt_chart::__tests__::test_support::tests::{internal_node, tip};
  use crate::rtt_chart::data::{RttChartData, RttPoint};
  use pretty_assertions::assert_eq;
  use treetime_utils::assert_error;

  #[test]
  fn test_rtt_chart_data_shows_dated_tips_and_counts_undated_ones() {
    let results = [
      tip(Some(2001.0), 0.5, false),
      internal_node(Some(2000.5), 0.25),
      tip(None, 0.75, false),
      tip(Some(2004.0), 2.0, true),
      tip(None, 1.0, true),
    ];
    let expected = RttChartData {
      fitted: vec![RttPoint { date: 2001.0, div: 0.5 }],
      excluded: vec![RttPoint { date: 2004.0, div: 2.0 }],
      undated: 2,
    };
    assert_eq!(expected, RttChartData::new(&results).unwrap());
  }

  #[test]
  fn test_rtt_chart_data_without_dated_tips_is_an_error() {
    assert_error!(
      RttChartData::new(&[tip(None, 0.5, false), internal_node(Some(2000.0), 0.0)]),
      "the root-to-tip chart needs at least one sample with a date"
    );
  }

  #[test]
  fn test_rtt_chart_data_fit_range_covers_the_fitted_tips() {
    let data = RttChartData::new(&[
      tip(Some(2003.0), 1.5, false),
      tip(Some(2001.0), 0.5, false),
      tip(Some(2009.0), 0.5, true),
    ])
    .unwrap();
    assert_eq!(
      ((2001.0, 2003.0), (2001.0, 2009.0)),
      (data.fit_range(), data.date_range())
    );
  }

  #[test]
  fn test_rtt_chart_data_fit_range_without_fitted_tips_covers_the_excluded_ones() {
    let data = RttChartData::new(&[tip(Some(2003.0), 1.5, true), tip(Some(2001.0), 0.5, true)]).unwrap();
    assert_eq!((2001.0, 2003.0), data.fit_range());
  }
}
