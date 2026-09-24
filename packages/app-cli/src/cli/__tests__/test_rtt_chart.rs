#[cfg(test)]
mod tests {
  use crate::cli::rtt_chart::write_clock_regression_chart_text;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime::clock::clock_model::{ClockModel, ClockRegression};
  use treetime::clock::rtt::ClockRegressionResult;
  use treetime_utils::io::json::json_read_str;

  #[test]
  fn test_rtt_chart_text_writes_table_then_chart() {
    let clock_model = helpers::clock_model();
    let results = vec![
      helpers::result("A", 2001.0, 0.002, false),
      helpers::result("B", 2003.0, 0.006, false),
      helpers::result("C", 2005.0, 0.030, true),
    ];

    let mut buf = Vec::new();
    write_clock_regression_chart_text(&mut buf, &results, &clock_model, (60, 20)).unwrap();
    let text = String::from_utf8(buf).unwrap();

    let (table, chart) = text.split_at(text.find('╯').unwrap() + '╯'.len_utf8() + 1);
    let x_axis_labels = chart.lines().rfind(|line| !line.is_empty()).unwrap();

    assert_eq!(
      indoc! {"
        ╭──────────────────┬──────────────────╮
        │ Clock regression │ div = 0.002t - 4 │
        ├──────────────────┼──────────────────┤
        │ tMRCA            │ 2000.0           │
        ├──────────────────┼──────────────────┤
        │ Rate             │ 0.002000         │
        ├──────────────────┼──────────────────┤
        │ Intercept        │ -4.0000          │
        ├──────────────────┼──────────────────┤
        │ R                │ 0.9900           │
        ├──────────────────┼──────────────────┤
        │ R²               │ 0.9801           │
        ├──────────────────┼──────────────────┤
        │ χ²               │ 0.000e0          │
        ╰──────────────────┴──────────────────╯
      "},
      table
    );
    assert_eq!("2001.0                  2005.0", x_axis_labels);
  }

  mod helpers {
    use super::*;

    pub(super) fn clock_model() -> ClockModel {
      let regression: ClockRegression = json_read_str(indoc! {r#"{
        "clock_rate": 0.002,
        "intercept": -4.0,
        "chisq": 0.0,
        "r_val": 0.99,
        "hessian": [[0.0, 0.0], [0.0, 0.0]],
        "cov": [[1e-8, 0.0], [0.0, 0.5]]
      }"#})
      .unwrap();
      ClockModel::from_regression(&regression).unwrap()
    }

    pub(super) fn result(name: &str, date: f64, div: f64, is_outlier: bool) -> ClockRegressionResult {
      ClockRegressionResult {
        name: Some(name.to_owned()),
        div,
        date: Some(date),
        predicted_date: date,
        clock_deviation: None,
        is_outlier,
        is_leaf: true,
      }
    }
  }
}
