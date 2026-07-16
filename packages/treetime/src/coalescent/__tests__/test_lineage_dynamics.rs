#[cfg(test)]
mod tests {
  use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
  use crate::coalescent::time_coordinate::CalendarTime;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;

  #[test]
  fn test_lineage_dynamics_binary_tree_calendar_direction() -> Result<(), Report> {
    // Oracle: one ancestral lineage exists before the root. Crossing a binary
    // merger toward the present creates two lineages; sampling removes both.
    let events = vec![(calendar(2000.0), -1), (calendar(2010.0), 1), (calendar(2010.0), 1)];

    let actual = compute_lineage_count_distribution(&events)?;

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

    let actual = compute_lineage_count_distribution(&events)?;

    pretty_assert_ulps_eq!(actual.eval_left(2000.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2000.0), 4.0, max_ulps = 4);
    pretty_assert_ulps_eq!(actual.eval(2010.0), 0.0, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_lineage_dynamics_rejects_incomplete_event_balance() {
    let events = vec![(calendar(2000.0), -1), (calendar(2010.0), 1)];

    let error = compute_lineage_count_distribution(&events).unwrap_err();

    assert!(error.to_string().contains("must end at zero"));
  }

  #[test]
  fn test_lineage_dynamics_rejects_empty_events() {
    let error = compute_lineage_count_distribution(&[]).unwrap_err();
    assert!(error.to_string().contains("empty events"));
  }

  fn calendar(time: f64) -> CalendarTime {
    CalendarTime::new(time)
  }
}
