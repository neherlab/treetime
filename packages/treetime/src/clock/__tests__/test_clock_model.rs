#[cfg(test)]
mod tests {
  use crate::clock::clock_model::{ClockLine, ClockModel, ClockModelStats, ClockRegression, RegressionStats};
  use crate::clock::clock_set::ClockSet;
  use eyre::Report;
  use indoc::indoc;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  fn clock_set_with_rate(target_rate: f64) -> ClockSet {
    let t1 = 0.0;
    let t2 = 10.0;
    let d1 = 0.0;
    let d2 = target_rate * (t2 - t1);
    let variance = 1.0;
    let cs1 = ClockSet::leaf_contribution_to_parent(Some(t1), d1, variance);
    let cs2 = ClockSet::leaf_contribution_to_parent(Some(t2), d2, variance);
    &cs1 + &cs2
  }

  #[test]
  fn test_clock_regression_positive_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::try_from(&cs)?;
    assert!(reg.clock_rate() > 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_regression_allows_negative_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(-0.005);
    let reg = ClockRegression::try_from(&cs)?;
    assert!(reg.clock_rate() < 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_regression_allows_zero_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::try_from(&cs)?;
    assert!((reg.clock_rate()).abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_regression_clock_deviation() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::try_from(&cs)?;
    let dev = reg.clock_deviation(2020.0, 0.5);
    let expected = 2020.0 * reg.clock_rate() + reg.intercept() - 0.5;
    assert!((dev - expected).abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_positive_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression(&reg)?;
    assert!(model.clock_rate() > 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_rejects_negative() {
    let cs = clock_set_with_rate(-0.005);
    let reg = ClockRegression::try_from(&cs).unwrap();
    let err = ClockModel::from_regression(&reg).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("non-positive"), "expected 'non-positive' in: {msg}");
    assert!(msg.contains("--clock-rate"), "expected '--clock-rate' in: {msg}");
  }

  #[test]
  fn test_clock_model_from_regression_rejects_zero() {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::try_from(&cs).unwrap();
    let err = ClockModel::from_regression(&reg).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("non-positive"), "expected 'non-positive' in: {msg}");
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_negative() -> Result<(), Report> {
    let cs = clock_set_with_rate(-0.005);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg);
    assert!((model.clock_rate() - reg.clock_rate()).abs() < 1e-15);
    assert!(model.clock_rate() < 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_zero() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg);
    assert!(model.clock_rate().abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_positive() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg);
    assert!((model.clock_rate() - reg.clock_rate()).abs() < 1e-15);
    assert!(model.clock_rate() > 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_model_with_fixed_rate_positive() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let model = ClockModel::with_fixed_rate(&cs, 0.001)?;
    assert!((model.clock_rate() - 0.001).abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_model_with_fixed_rate_rejects_negative() {
    let cs = clock_set_with_rate(0.003);
    let err = ClockModel::with_fixed_rate(&cs, -0.001).unwrap_err();
    let msg = err.to_string();
    assert!(
      msg.contains("must be positive"),
      "expected 'must be positive' in: {msg}"
    );
  }

  #[test]
  fn test_clock_model_with_fixed_rate_rejects_zero() {
    let cs = clock_set_with_rate(0.003);
    let err = ClockModel::with_fixed_rate(&cs, 0.0).unwrap_err();
    let msg = err.to_string();
    assert!(
      msg.contains("must be positive"),
      "expected 'must be positive' in: {msg}"
    );
  }

  #[test]
  fn test_clock_line_deviation_consistent_between_types() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression(&reg)?;
    let date = 2020.0;
    let div = 0.5;
    assert!((reg.clock_deviation(date, div) - model.clock_deviation(date, div)).abs() < 1e-15);
    Ok(())
  }

  fn estimated_model() -> ClockModel {
    let stats = ClockModelStats::Estimated(RegressionStats {
      chisq: 1.5,
      r_val: 0.9,
      hessian: array![[1.0, 2.0], [3.0, 4.0]],
      cov: array![[0.1, 0.0], [0.0, 0.2]],
    });
    ClockModel::for_testing_with_stats(0.003, -6.0, stats)
  }

  #[test]
  fn test_clock_model_serializes_stats_matrices_as_nested_arrays() -> Result<(), Report> {
    let expected = indoc! {r#"{
      "clock_rate": 0.003,
      "intercept": -6.0,
      "stats": {
        "estimated": {
          "chisq": 1.5,
          "r_val": 0.9,
          "hessian": [
            [
              1.0,
              2.0
            ],
            [
              3.0,
              4.0
            ]
          ],
          "cov": [
            [
              0.1,
              0.0
            ],
            [
              0.0,
              0.2
            ]
          ]
        }
      }
    }"#};
    let actual = json_write_str(&estimated_model(), JsonPretty(true))?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_clock_model_stats_matrices_json_roundtrip() -> Result<(), Report> {
    let json = json_write_str(&estimated_model(), JsonPretty(true))?;
    let restored: ClockModel = serde_json::from_str(&json)?;
    let ClockModelStats::Estimated(stats) = restored.stats() else {
      panic!("expected Estimated stats after round-trip");
    };
    assert_eq!(&array![[1.0, 2.0], [3.0, 4.0]], &stats.hessian);
    assert_eq!(&array![[0.1, 0.0], [0.0, 0.2]], &stats.cov);
    Ok(())
  }
}
