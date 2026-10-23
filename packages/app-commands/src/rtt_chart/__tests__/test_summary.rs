#[cfg(test)]
mod tests {
  use crate::rtt_chart::__tests__::test_support::tests::{
    estimated_model, fitted_tips, fixed_model, model_with_rate_std, tip,
  };
  use crate::rtt_chart::data::RttChartData;
  use crate::rtt_chart::summary::RttSummary;
  use helpers::statistic;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::iter::repeat_with;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_rtt_summary_states_rate_root_date_and_fit_quality() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.002), (2005.0, 0.01)])).unwrap();
    let expected = RttSummary {
      rate: statistic("Clock rate", "2.00 × 10⁻³", Some("substitutions per site per year")),
      root_date: statistic("Root date", "2000.00", Some("January 2000")),
      r_squared: Some(statistic("R²", "0.980", None)),
      samples: "2 samples in the fit".to_owned(),
      notes: vec![],
    };
    assert_eq!(expected, RttSummary::new(&data, &estimated_model(0.002, -4.0, 0.99)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::months_resolve_short_spans(  0.2,   ("2000.000", Some("January 2000")))]
  #[case::two_decimals_below_ten_years(4.0,   ("2000.00",  Some("January 2000")))]
  #[case::one_decimal_for_decades(     40.0,  ("2000.0",   None))]
  #[case::whole_years_for_centuries(   200.0, ("2000",     None))]
  #[trace]
  fn test_rtt_summary_root_date_precision_follows_the_sampling_span(
    #[case] span: f64,
    #[case] (value, detail): (&str, Option<&str>),
  ) {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.002), (2001.0 + span, 0.01)])).unwrap();
    let summary = RttSummary::new(&data, &estimated_model(0.002, -4.0, 0.99));
    assert_eq!(statistic("Root date", value, detail), summary.root_date);
  }

  #[test]
  fn test_rtt_summary_root_date_is_where_the_fitted_line_reaches_zero_divergence() {
    let data = RttChartData::new(&fitted_tips(&[(2017.77, 4.7e-5), (2024.62, 4.7e-4)])).unwrap();
    let summary = RttSummary::new(&data, &estimated_model(5.072_949e-5, -0.102_306_463_7, 0.96));
    assert_eq!("2016.71", summary.root_date.value);
  }

  #[test]
  fn test_rtt_summary_without_positive_rate_has_no_root_date() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.01), (2005.0, 0.002)])).unwrap();
    let summary = RttSummary::new(&data, &estimated_model(-0.002, 4.012, -0.5));
    assert_eq!(
      (
        statistic(
          "Clock rate",
          "\u{2212}2.00 × 10⁻³",
          Some("substitutions per site per year")
        ),
        statistic("Root date", "None", Some("no positive clock rate")),
      ),
      (summary.rate, summary.root_date)
    );
  }

  #[test]
  fn test_rtt_summary_of_a_fixed_rate_has_no_fit_quality() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.002), (2005.0, 0.01)])).unwrap();
    let summary = RttSummary::new(&data, &fixed_model(0.002, -4.0));
    assert_eq!(
      ("Clock rate (fixed)", None),
      (summary.rate.label.as_str(), summary.r_squared)
    );
  }

  #[test]
  fn test_rtt_summary_rate_carries_its_standard_deviation() {
    let data = RttChartData::new(&fitted_tips(&[(2001.0, 0.002), (2005.0, 0.01)])).unwrap();
    let summary = RttSummary::new(&data, &model_with_rate_std(0.002, -4.0, 1.5e-4));
    assert_eq!("2.00 ± 0.15 × 10⁻³", summary.rate.value);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::all_fitted(       (3, 0, 0), ("3 samples in the fit",                                vec![]))]
  #[case::one_fitted(       (1, 0, 0), ("1 sample in the fit",                                 vec![]))]
  #[case::with_excluded(    (3, 2, 0), ("3 samples in the fit, 2 excluded by the clock filter", vec![]))]
  #[case::one_undated(      (3, 0, 1), ("3 samples in the fit",                                vec!["1 sample without a date is not shown"]))]
  #[case::several_undated(  (3, 1, 4), ("3 samples in the fit, 1 excluded by the clock filter", vec!["4 samples without a date are not shown"]))]
  #[trace]
  fn test_rtt_summary_counts_the_samples(
    #[case] (fitted, excluded, undated): (usize, usize, usize),
    #[case] (samples, notes): (&str, Vec<&str>),
  ) {
    let results = (0..fitted)
      .map(|index| tip(Some(2000.0 + index.to_f64_exact()), 0.01, false))
      .chain(repeat_with(|| tip(Some(2003.0), 0.05, true)).take(excluded))
      .chain(repeat_with(|| tip(None, 0.01, false)).take(undated))
      .collect::<Vec<_>>();
    let summary = RttSummary::new(&RttChartData::new(&results).unwrap(), &estimated_model(0.002, -4.0, 0.9));
    assert_eq!((samples, notes), (summary.samples.as_str(), summary.notes.iter().map(String::as_str).collect()));
  }

  #[test]
  fn test_rtt_summary_description_reads_as_sentences() {
    let data = RttChartData::new(&[
      tip(Some(2001.0), 0.002, false),
      tip(Some(2005.0), 0.01, false),
      tip(None, 0.01, false),
    ])
    .unwrap();
    let summary = RttSummary::new(&data, &estimated_model(0.002, -4.0, 0.99));
    assert_eq!(
      "Clock rate 2.00 × 10⁻³ (substitutions per site per year). Root date 2000.00 (January 2000). R² 0.980. \
       2 samples in the fit. 1 sample without a date is not shown",
      summary.description()
    );
  }

  mod helpers {
    use crate::rtt_chart::summary::Statistic;

    pub(super) fn statistic(label: &str, value: &str, detail: Option<&str>) -> Statistic {
      Statistic {
        label: label.to_owned(),
        value: value.to_owned(),
        detail: detail.map(str::to_owned),
      }
    }
  }
}
