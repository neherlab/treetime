#[cfg(test)]
pub(crate) mod tests {
  use indoc::formatdoc;
  use treetime::clock::clock_model::ClockModel;
  use treetime::clock::rtt::ClockRegressionResult;
  use treetime_utils::io::json::json_read_str;

  pub(crate) fn tip(date: Option<f64>, div: f64, is_outlier: bool) -> ClockRegressionResult {
    ClockRegressionResult {
      name: String::new(),
      div,
      date,
      predicted_date: 0.0,
      clock_deviation: None,
      is_outlier,
      is_leaf: true,
      date_source: None,
    }
  }

  pub(crate) fn internal_node(date: Option<f64>, div: f64) -> ClockRegressionResult {
    ClockRegressionResult {
      is_leaf: false,
      ..tip(date, div, false)
    }
  }

  pub(crate) fn fitted_tips(points: &[(f64, f64)]) -> Vec<ClockRegressionResult> {
    points.iter().map(|&(date, div)| tip(Some(date), div, false)).collect()
  }

  pub(crate) fn estimated_model(clock_rate: f64, intercept: f64, r_val: f64) -> ClockModel {
    json_read_str(formatdoc! {r#"{{
      "clock_rate": {clock_rate:e},
      "intercept": {intercept:e},
      "stats": {{ "estimated": {{ "chisq": 0.0, "r_val": {r_val:e} }} }}
    }}"#})
    .unwrap()
  }

  pub(crate) fn model_with_rate_std(clock_rate: f64, intercept: f64, rate_std: f64) -> ClockModel {
    json_read_str(formatdoc! {r#"{{
      "clock_rate": {clock_rate:e},
      "intercept": {intercept:e},
      "stats": {{ "estimated": {{ "chisq": 0.0, "r_val": 0.5, "rate_std": {rate_std:e} }} }}
    }}"#})
    .unwrap()
  }

  pub(crate) fn fixed_model(clock_rate: f64, intercept: f64) -> ClockModel {
    json_read_str(formatdoc! {r#"{{
      "clock_rate": {clock_rate:e},
      "intercept": {intercept:e},
      "stats": "fixed"
    }}"#})
    .unwrap()
  }
}
