#[cfg(test)]
mod tests {
  use crate::timetree::convergence::metrics::{ConvergenceMetrics, NODE_TIME_TOLERANCE_YEARS};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_primitives::LogLh;
  use treetime_utils::io::json::{JsonPretty, json_write_str};

  #[test]
  fn test_metrics_serializes_explicit_log_lh_schema() -> Result<(), Report> {
    let metrics = ConvergenceMetrics {
      n_diff: 2,
      n_resolved: 1,
      max_time_change: Some(0.5),
      rms_time_change: Some(0.125),
      log_lh_seq: Some(LogLh::new(-10.0)),
      log_lh_time: Some(LogLh::new(-20.0)),
      log_lh_coal: Some(LogLh::new(-30.0)),
    };
    let expected = indoc! {r#"{
      "n_diff": 2,
      "n_resolved": 1,
      "max_time_change": 0.5,
      "rms_time_change": 0.125,
      "log_lh_seq": -10.0,
      "log_lh_time": -20.0,
      "log_lh_coal": -30.0
    }"#};

    let actual = json_write_str(&metrics, JsonPretty(true))?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_time_change_and_no_differences(      (0, 0, None,                                   None),                                   true)]
  #[case::no_time_change_and_differences(         (2, 0, None,                                   None),                                   false)]
  #[case::settled_times_override_differences(     (2, 0, Some(0.1 * NODE_TIME_TOLERANCE_YEARS),  Some(0.1 * NODE_TIME_TOLERANCE_YEARS)),  true)]
  #[case::times_move_by_the_tolerance(            (0, 0, Some(NODE_TIME_TOLERANCE_YEARS),        Some(NODE_TIME_TOLERANCE_YEARS)),        false)]
  #[case::times_move_beyond_the_tolerance(        (0, 0, Some(10.0 * NODE_TIME_TOLERANCE_YEARS), Some(10.0 * NODE_TIME_TOLERANCE_YEARS)), false)]
  #[case::max_moves_while_rms_settles(            (0, 0, Some(NODE_TIME_TOLERANCE_YEARS),        Some(0.1 * NODE_TIME_TOLERANCE_YEARS)),  false)]
  #[case::polytomies_resolve_with_settled_times(  (0, 1, Some(0.1 * NODE_TIME_TOLERANCE_YEARS),  Some(0.1 * NODE_TIME_TOLERANCE_YEARS)),  false)]
  #[case::polytomies_resolve_without_time_change( (0, 1, None,                                   None),                                   false)]
  #[trace]
  fn test_metrics_convergence_requires_settled_times_and_no_resolved_polytomies(
    #[case] (n_diff, n_resolved, max_time_change, rms_time_change): (usize, usize, Option<f64>, Option<f64>),
    #[case] expected: bool,
  ) {
    let metrics = ConvergenceMetrics {
      n_diff,
      n_resolved,
      max_time_change,
      rms_time_change,
      log_lh_seq: None,
      log_lh_time: None,
      log_lh_coal: None,
    };

    assert_eq!(expected, metrics.has_converged());
  }
}
