#[cfg(test)]
mod tests {
  use crate::coalescent::lineage_dynamics::compute_lineage_count_distribution;
  use crate::coalescent::time_coordinate::Tbp;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;

  #[test]
  fn test_piecewise_constant_eval() {
    // Breakpoints: [1.0, 5.0, 10.0]
    // Values: [0.0, 1.0, 2.0, 3.0]
    // t < 1.0 -> 0.0
    // 1.0 <= t < 5.0 -> 1.0
    // 5.0 <= t < 10.0 -> 2.0
    // t >= 10.0 -> 3.0
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_ulps_eq!(pc.eval(-1.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(0.5), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(1.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(3.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(5.0), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(7.5), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(10.0), 3.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(100.0), 3.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_constant_eval_many() {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0], array![0.0, 1.0, 2.0]);
    let ts = array![0.0, 1.0, 3.0, 5.0, 10.0];
    let result = pc.eval_many(&ts);
    pretty_assert_ulps_eq!(result[0], 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[1], 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[2], 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[3], 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[4], 2.0, max_ulps = 4);
  }

  fn tbp(t: f64) -> Tbp {
    Tbp::new(t)
  }

  #[test]
  fn test_lineage_count_simple_tree() -> Result<(), Report> {
    // Simple tree: 2 tips at 0, root at 10.
    // Events:
    // 0.0: Tip 1 (+1)
    // 0.0: Tip 2 (+1)
    // 10.0: Root (-1) (merger of 2 -> 1)
    let events = vec![(tbp(0.0), 1), (tbp(0.0), 1), (tbp(10.0), -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Check values at specific points
    // Before any events (t < 0.0): 0
    pretty_assert_ulps_eq!(lineage_counts.eval(-1.0), 0.0, max_ulps = 4);
    // After t=0 events (+2): 2
    pretty_assert_ulps_eq!(lineage_counts.eval(5.0), 2.0, max_ulps = 4);
    // After t=10 event (-1): 1
    pretty_assert_ulps_eq!(lineage_counts.eval(15.0), 1.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_lineage_count_single_event() -> Result<(), Report> {
    let events = vec![(tbp(5.0), 2)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Before event (t < 5.0): 0
    pretty_assert_ulps_eq!(lineage_counts.eval(4.0), 0.0, max_ulps = 4);
    // After event (t >= 5.0): 2
    pretty_assert_ulps_eq!(lineage_counts.eval(6.0), 2.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_lineage_count_empty_events() {
    let events: Vec<(Tbp, i32)> = vec![];
    drop(compute_lineage_count_distribution(&events).unwrap_err());
  }

  #[test]
  fn test_lineage_count_aggregation() -> Result<(), Report> {
    let events = vec![(tbp(5.0), 1), (tbp(5.0), 1), (tbp(5.0), -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Net delta +1 at 5.0
    // Before 5.0: 0
    pretty_assert_ulps_eq!(lineage_counts.eval(4.5), 0.0, max_ulps = 4);
    // After 5.0: 1
    pretty_assert_ulps_eq!(lineage_counts.eval(5.5), 1.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_lineage_count_decreasing_then_increasing() -> Result<(), Report> {
    let events = vec![(tbp(0.0), 3), (tbp(5.0), -1), (tbp(10.0), 1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // 0.0: +3 -> 3
    // 5.0: -1 -> 2
    // 10.0: +1 -> 3
    // Before first event: 0
    pretty_assert_ulps_eq!(lineage_counts.eval(-1.0), 0.0, max_ulps = 4);
    // After 0.0, before 5.0: 3
    pretty_assert_ulps_eq!(lineage_counts.eval(2.5), 3.0, max_ulps = 4);
    // After 5.0, before 10.0: 2
    pretty_assert_ulps_eq!(lineage_counts.eval(7.5), 2.0, max_ulps = 4);
    // After 10.0: 3
    pretty_assert_ulps_eq!(lineage_counts.eval(15.0), 3.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_lineage_count_breakpoints() -> Result<(), Report> {
    let events = vec![(tbp(10.0), 1), (tbp(20.0), -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let breakpoints = lineage_counts.breakpoints();
    assert_eq!(breakpoints.len(), 2);
    pretty_assert_ulps_eq!(breakpoints[0], 10.0, max_ulps = 4);
    pretty_assert_ulps_eq!(breakpoints[1], 20.0, max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_lineage_count_negative_deltas() -> Result<(), Report> {
    let events = vec![(tbp(0.0), 5), (tbp(5.0), -1), (tbp(10.0), -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // 0.0: +5 -> 5
    // 5.0: -1 -> 4
    // 10.0: -1 -> 3
    // Before first event: 0
    pretty_assert_ulps_eq!(lineage_counts.eval(-1.0), 0.0, max_ulps = 4);
    // After 0.0, before 5.0: 5
    pretty_assert_ulps_eq!(lineage_counts.eval(2.5), 5.0, max_ulps = 4);
    // After 5.0, before 10.0: 4
    pretty_assert_ulps_eq!(lineage_counts.eval(7.5), 4.0, max_ulps = 4);
    // After 10.0: 3
    pretty_assert_ulps_eq!(lineage_counts.eval(15.0), 3.0, max_ulps = 4);

    Ok(())
  }
}
