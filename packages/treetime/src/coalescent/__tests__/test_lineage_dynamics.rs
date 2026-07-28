#[cfg(test)]
mod tests {
  use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
  use crate::coalescent::time_coordinate::CalendarTime;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use treetime_utils::assert_error;

  #[test]
  fn test_lineage_dynamics_binary_tree_calendar_direction() -> Result<(), Report> {
    // Oracle: one ancestral lineage exists before the root. Crossing a binary
    // merger toward the present creates two lineages; sampling removes both.
    let events = vec![(calendar(2000.0), -1), (calendar(2010.0), 1), (calendar(2010.0), 1)];

    let actual = compute_lineage_count_distribution(&events, 0)?;

    pretty_assert_ulps_eq!(actual.eval(1999.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2000.0), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2005.0), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2010.0), 0.0, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_lineage_dynamics_polytomy_uses_child_count_minus_one() -> Result<(), Report> {
    // A four-child event changes the lineage count by three.
    let events = vec![
      (calendar(2000.0), -3),
      (calendar(2010.0), 1),
      (calendar(2010.0), 1),
      (calendar(2010.0), 1),
      (calendar(2010.0), 1),
    ];

    let actual = compute_lineage_count_distribution(&events, 0)?;

    pretty_assert_ulps_eq!(actual.eval_left(2000.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2000.0), 4.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2010.0), 0.0, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_lineage_dynamics_rejects_incomplete_event_balance() {
    let events = vec![(calendar(2000.0), -1), (calendar(2010.0), 1)];

    assert_error!(
      compute_lineage_count_distribution(&events, 0),
      "Lineage count must end at 0 after the latest retained sample, got 1"
    );
  }

  #[test]
  fn test_lineage_dynamics_rejects_empty_events() {
    assert_error!(
      compute_lineage_count_distribution(&[], 0),
      "Cannot build lineage count from empty events"
    );
  }

  #[test]
  fn test_lineage_dynamics_retains_lineage_for_excluded_subtree() -> Result<(), Report> {
    // Oracle: filtering one bad leaf from a three-child root leaves one lineage
    // after the latest retained sample in v0's `calc_branch_count()`.
    let events = vec![(calendar(2000.0), -2), (calendar(2005.0), 1), (calendar(2010.0), 1)];

    let actual = compute_lineage_count_distribution(&events, 1)?;

    pretty_assert_ulps_eq!(actual.eval(1999.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2000.0), 3.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2005.0), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2010.0), 1.0, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_lineage_dynamics_rejects_negative_terminal_count() {
    let events = vec![(calendar(2000.0), -1), (calendar(2010.0), 1)];

    assert_error!(
      compute_lineage_count_distribution(&events, -1),
      "Terminal lineage count must be non-negative, got -1"
    );
  }

  fn calendar(time: f64) -> CalendarTime {
    CalendarTime::new(time)
  }
}
