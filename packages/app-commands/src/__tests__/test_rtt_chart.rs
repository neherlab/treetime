#[cfg(test)]
mod tests {
  use crate::rtt_chart::gather_points;
  use helpers::{clock_model, row};
  use pretty_assertions::assert_eq;
  use treetime_utils::assert_error;

  #[test]
  fn test_rtt_chart_without_dated_samples_is_an_error() {
    assert_error!(
      gather_points(&[row(None, false), row(None, true)], &clock_model()),
      "the root-to-tip chart needs at least one sample with a date"
    );
  }

  #[test]
  fn test_rtt_chart_line_spans_the_dated_samples() {
    let points = gather_points(
      &[row(Some(2000.0), false), row(Some(2010.0), true), row(None, false)],
      &clock_model(),
    )
    .unwrap();
    assert_eq!(
      ((2000.0, 0.0), (2010.0, 0.02), 1, 1),
      (
        points.line[0],
        points.line[1],
        points.norm_points.len(),
        points.outlier_points.len()
      )
    );
  }

  mod helpers {
    use indoc::indoc;
    use treetime::clock::clock_model::ClockModel;
    use treetime::clock::rtt::ClockRegressionResult;
    use treetime_utils::io::json::json_read_str;

    pub(super) fn row(date: Option<f64>, is_outlier: bool) -> ClockRegressionResult {
      ClockRegressionResult {
        name: None,
        div: 0.01,
        date,
        predicted_date: 2005.0,
        clock_deviation: None,
        is_outlier,
        is_leaf: true,
        date_source: None,
      }
    }

    pub(super) fn clock_model() -> ClockModel {
      json_read_str(indoc! {r#"{
        "clock_rate": 0.002,
        "intercept": -4.0,
        "stats": {
          "estimated": {
            "chisq": 0.0,
            "r_val": 0.99,
            "hessian": [[0.0, 0.0], [0.0, 0.0]],
            "cov": [[1e-8, 0.0], [0.0, 0.5]]
          }
        }
      }"#})
      .unwrap()
    }
  }
}
