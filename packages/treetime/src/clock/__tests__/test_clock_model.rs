#[cfg(test)]
mod tests {
  use crate::clock::clock_model::{ClockLine, ClockModel, ClockModelStats, ClockRegression, RegressionStats};
  use crate::clock::clock_set::ClockSet;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use indoc::indoc;

  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::json_read_str;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  fn clock_set_with_rate(target_rate: f64) -> ClockSet {
    let t1 = 0.0;
    let t2 = 10.0;
    let d1 = 0.0;
    let d2 = target_rate * (t2 - t1);
    let cs1 = ClockSet::leaf_contribution(Some(t1)).shift(d1);
    let cs2 = ClockSet::leaf_contribution(Some(t2)).shift(d2);
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
    let model = ClockModel::from_regression_allow_negative(&reg, &NoopProgress);
    assert!((model.clock_rate() - reg.clock_rate()).abs() < 1e-15);
    assert!(model.clock_rate() < 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_zero() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg, &NoopProgress);
    assert!(model.clock_rate().abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_positive() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::try_from(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg, &NoopProgress);
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

  fn estimated_model(rate_std: Option<f64>) -> ClockModel {
    let stats = ClockModelStats::Estimated(RegressionStats {
      chisq: 1.5,
      r_val: 0.9,
      rate_std,
    });
    ClockModel::for_testing_with_stats(0.003, -6.0, stats)
  }

  #[test]
  fn test_clock_model_serializes_the_rate_std_only_when_known() -> Result<(), Report> {
    let with_std = indoc! {r#"{
      "clock_rate": 0.003,
      "intercept": -6.0,
      "stats": {
        "estimated": {
          "chisq": 1.5,
          "r_val": 0.9,
          "rate_std": 0.0002
        }
      }
    }"#};
    let without_std = indoc! {r#"{
      "clock_rate": 0.003,
      "intercept": -6.0,
      "stats": {
        "estimated": {
          "chisq": 1.5,
          "r_val": 0.9
        }
      }
    }"#};
    assert_eq!(
      (with_std, without_std),
      (
        json_write_str(&estimated_model(Some(2e-4)), JsonPretty(true))?.as_str(),
        json_write_str(&estimated_model(None), JsonPretty(true))?.as_str()
      )
    );
    Ok(())
  }

  #[test]
  fn test_clock_model_stats_json_roundtrip() -> Result<(), Report> {
    for rate_std in [Some(2e-4), None] {
      let json = json_write_str(&estimated_model(rate_std), JsonPretty(true))?;
      let restored: ClockModel = json_read_str(&json)?;
      assert_eq!(
        (rate_std, Some(1.5), Some(0.9)),
        (restored.rate_std(), restored.chisq(), restored.r_val())
      );
    }
    Ok(())
  }

  #[test]
  fn test_clock_model_estimated_at_the_regression_rate_reproduces_the_regression_stats() -> Result<(), Report> {
    let set = ClockSet::leaf_contribution(Some(2019.0)).shift(0.01)
      + ClockSet::leaf_contribution(Some(2020.5)).shift(0.03)
      + ClockSet::leaf_contribution(Some(2021.0)).shift(0.05);
    let regression = ClockRegression::try_from(&set)?;
    let model = ClockModel::from_regression(&regression)?;

    let restated = model.estimated_at_rate(&set, Some(1e-3));

    approx::assert_relative_eq!(model.chisq().unwrap(), restated.chisq().unwrap(), max_relative = 1e-6);
    assert_eq!((model.r_val(), Some(1e-3)), (restated.r_val(), restated.rate_std()));
    Ok(())
  }
}
