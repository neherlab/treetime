#[cfg(test)]
mod tests {
  use crate::timetree::convergence::metrics::ConvergenceMetrics;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  #[test]
  fn test_metrics_serializes_explicit_log_lh_schema() -> Result<(), Report> {
    let metrics = ConvergenceMetrics {
      n_diff: 2,
      n_resolved: 1,
      log_lh_seq: Some(-10.0),
      log_lh_pos: Some(-20.0),
      log_lh_coal: Some(-30.0),
      log_lh_total: Some(-60.0),
    };
    let expected = indoc! {r#"{
      "n_diff": 2,
      "n_resolved": 1,
      "log_lh_seq": -10.0,
      "log_lh_pos": -20.0,
      "log_lh_coal": -30.0,
      "log_lh_total": -60.0
    }"#};

    let actual = json_write_str(&metrics, JsonPretty(true))?;

    assert_eq!(expected, actual);
    Ok(())
  }
}
