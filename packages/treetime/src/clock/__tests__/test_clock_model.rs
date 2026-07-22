#[cfg(test)]
mod tests {
  use crate::clock::clock_model::{ClockLine, ClockModel, ClockRegression};
  use crate::payload::clock_set::ClockSet;
  use eyre::Report;

  /// Build a ClockSet from two leaves whose regression produces the given clock_rate.
  ///
  /// Uses leaf_contribution_to_parent with unit variance, placing leaves at
  /// dates t1=0, t2=10 with divergences chosen to produce the target rate.
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

  // --- ClockRegression ---

  #[test]
  fn test_clock_regression_positive_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::from_clock_set(&cs)?;
    assert!(reg.clock_rate() > 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_regression_allows_negative_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(-0.005);
    let reg = ClockRegression::from_clock_set(&cs)?;
    assert!(reg.clock_rate() < 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_regression_allows_zero_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::from_clock_set(&cs)?;
    assert!((reg.clock_rate()).abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_regression_clock_deviation() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::from_clock_set(&cs)?;
    let dev = reg.clock_deviation(2020.0, 0.5);
    let expected = 2020.0 * reg.clock_rate() + reg.intercept() - 0.5;
    assert!((dev - expected).abs() < 1e-15);
    Ok(())
  }

  // --- ClockModel::from_regression ---

  #[test]
  fn test_clock_model_from_regression_positive_rate() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::from_clock_set(&cs)?;
    let model = ClockModel::from_regression(&reg)?;
    assert!(model.clock_rate() > 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_rejects_negative() {
    let cs = clock_set_with_rate(-0.005);
    let reg = ClockRegression::from_clock_set(&cs).unwrap();
    let err = ClockModel::from_regression(&reg).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("non-positive"), "expected 'non-positive' in: {msg}");
    assert!(msg.contains("--clock-rate"), "expected '--clock-rate' in: {msg}");
  }

  #[test]
  fn test_clock_model_from_regression_rejects_zero() {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::from_clock_set(&cs).unwrap();
    let err = ClockModel::from_regression(&reg).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("non-positive"), "expected 'non-positive' in: {msg}");
  }

  // --- ClockModel::from_regression_allow_negative ---
  // The clock command reports the regression without time inference, so a non-positive
  // rate is a valid result: the lenient constructor builds a model (warning) instead of
  // erroring. See kb/decisions/timetree-rejects-negative-clock-rate.md.

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_negative() -> Result<(), Report> {
    let cs = clock_set_with_rate(-0.005);
    let reg = ClockRegression::from_clock_set(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg);
    assert!((model.clock_rate() - reg.clock_rate()).abs() < 1e-15);
    assert!(model.clock_rate() < 0.0);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_zero() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.0);
    let reg = ClockRegression::from_clock_set(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg);
    assert!(model.clock_rate().abs() < 1e-15);
    Ok(())
  }

  #[test]
  fn test_clock_model_from_regression_allow_negative_builds_positive() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::from_clock_set(&cs)?;
    let model = ClockModel::from_regression_allow_negative(&reg);
    assert!((model.clock_rate() - reg.clock_rate()).abs() < 1e-15);
    assert!(model.clock_rate() > 0.0);
    Ok(())
  }

  // --- ClockModel::with_fixed_rate ---

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

  // --- ClockLine trait ---

  #[test]
  fn test_clock_line_deviation_consistent_between_types() -> Result<(), Report> {
    let cs = clock_set_with_rate(0.003);
    let reg = ClockRegression::from_clock_set(&cs)?;
    let model = ClockModel::from_regression(&reg)?;
    let date = 2020.0;
    let div = 0.5;
    assert!((reg.clock_deviation(date, div) - model.clock_deviation(date, div)).abs() < 1e-15);
    Ok(())
  }
}
