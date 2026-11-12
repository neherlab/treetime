#[cfg(test)]
mod tests {
  use crate::results::regression::{RegressionStatisticKind, StatisticValue, regression_statistics};
  use crate::results::year_date::YearDate;
  use helpers::{estimated_model, fixed_model, point, points, values};
  use pretty_assertions::assert_eq;
  use treetime_utils::pretty_assert_relative_eq;

  #[test]
  fn test_regression_statistics_match_ordinary_least_squares_of_the_regression_samples() {
    let mean_date = 2001.5;
    let mean_div = 0.002;
    let sxx: f64 = 5.0;
    let sxy = 0.006;
    let syy = 8e-6;
    let rate = sxy / sxx;
    let intercept = mean_div - rate * mean_date;
    let r = sxy / (sxx * syy).sqrt();
    let residual_sum_of_squares = syy - sxy * sxy / sxx;
    let model = estimated_model(rate, intercept, r);

    let statistics = regression_statistics(&model, &points());

    assert_eq!(
      vec![
        RegressionStatisticKind::ClockRate,
        RegressionStatisticKind::RootDate,
        RegressionStatisticKind::CorrelationCoefficient,
        RegressionStatisticKind::RSquared,
        RegressionStatisticKind::ResidualMeanSquare,
      ],
      statistics
        .iter()
        .map(|statistic| statistic.statistic)
        .collect::<Vec<_>>()
    );
    let expected = [
      0.0012,
      mean_date - mean_div / rate,
      3.0 / 10.0_f64.sqrt(),
      0.9,
      residual_sum_of_squares / 2.0,
    ];
    let actual = values(&statistics);
    pretty_assert_relative_eq!(expected.as_slice(), actual.as_slice(), max_relative = 1e-9);
  }

  #[test]
  fn test_regression_statistics_root_date_is_a_calendar_date() {
    let model = estimated_model(1e-3, -2.0, 0.9);

    let statistics = regression_statistics(&model, &points());

    assert_eq!(
      StatisticValue::Date {
        date: YearDate::new(2000.0)
      },
      statistics[1].value
    );
  }

  #[test]
  fn test_regression_statistics_of_a_fixed_rate_fit_one_parameter_and_omit_the_correlation() {
    let model = fixed_model(1e-3, -2.0);
    let points = [
      point(Some(2000.0), 0.001, false),
      point(Some(2001.0), 0.0005, false),
      point(Some(2002.0), 0.002, false),
    ];
    let residual_sum_of_squares = 1e-6 + 0.25e-6 + 0.0;

    let statistics = regression_statistics(&model, &points);

    assert_eq!(
      vec![
        RegressionStatisticKind::ClockRate,
        RegressionStatisticKind::RootDate,
        RegressionStatisticKind::ResidualMeanSquare,
      ],
      statistics
        .iter()
        .map(|statistic| statistic.statistic)
        .collect::<Vec<_>>()
    );
    let expected = [1e-3, 2000.0, residual_sum_of_squares / 2.0];
    let actual = values(&statistics);
    pretty_assert_relative_eq!(expected.as_slice(), actual.as_slice(), max_relative = 1e-9);
  }

  #[test]
  fn test_regression_statistics_omit_the_residual_mean_square_without_residual_degrees_of_freedom() {
    let model = estimated_model(1e-3, -2.0, 1.0);
    let points = [point(Some(2000.0), 0.0, false), point(Some(2001.0), 0.001, false)];

    let statistics = regression_statistics(&model, &points);

    assert_eq!(
      vec![
        RegressionStatisticKind::ClockRate,
        RegressionStatisticKind::RootDate,
        RegressionStatisticKind::CorrelationCoefficient,
        RegressionStatisticKind::RSquared,
      ],
      statistics
        .iter()
        .map(|statistic| statistic.statistic)
        .collect::<Vec<_>>()
    );
  }

  mod helpers {
    use crate::results::clock::RootToTipPoint;
    use crate::results::regression::{RegressionStatistic, StatisticValue};
    use crate::results::year_date::YearDate;
    use treetime::clock::clock_model::ClockModel;
    use treetime::clock::rtt::ClockDateSource;
    use treetime::o;
    use treetime_utils::io::json::json_read_str;

    pub(super) fn estimated_model(rate: f64, intercept: f64, r: f64) -> ClockModel {
      json_read_str(format!(
        r#"{{ "clock_rate": {rate:e}, "intercept": {intercept:e}, "stats": {{ "estimated": {{ "chisq": 0.0, "r_val": {r:e} }} }} }}"#
      ))
      .unwrap()
    }

    pub(super) fn fixed_model(rate: f64, intercept: f64) -> ClockModel {
      json_read_str(format!(
        r#"{{ "clock_rate": {rate:e}, "intercept": {intercept:e}, "stats": "fixed" }}"#
      ))
      .unwrap()
    }

    pub(super) fn points() -> Vec<RootToTipPoint> {
      vec![
        point(Some(2000.0), 0.0, false),
        point(Some(2001.0), 0.002, false),
        point(Some(2002.0), 0.002, false),
        point(Some(2003.0), 0.004, false),
        point(Some(2002.5), 0.05, true),
        point(None, 0.01, false),
      ]
    }

    pub(super) fn point(date: Option<f64>, div: f64, outlier: bool) -> RootToTipPoint {
      RootToTipPoint {
        name: o!("sample"),
        date: date.map(YearDate::new),
        date_source: Some(date.map_or(ClockDateSource::Missing, |_| ClockDateSource::Input)),
        div,
        predicted_date: YearDate::new(2000.0),
        residual_days: None,
        outlier,
      }
    }

    pub(super) fn values(statistics: &[RegressionStatistic]) -> Vec<f64> {
      statistics
        .iter()
        .map(|statistic| match &statistic.value {
          StatisticValue::Number { value } => *value,
          StatisticValue::Date { date } => date.year,
        })
        .collect()
    }
  }
}
