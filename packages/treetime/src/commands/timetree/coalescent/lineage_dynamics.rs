use crate::commands::timetree::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_utils::make_error;

/// Computes k(t) distribution from tree events.
///
/// k(t) is the number of concurrent lineages at time t.
/// The function is piecewise constant, stepping at each merger event.
/// Events must be sorted by increasing time (past to present).
pub fn compute_lineage_count_distribution(events: &[(f64, i32)]) -> Result<PiecewiseConstantFn, Report> {
  if events.is_empty() {
    return make_error!("Cannot build lineage count from empty events");
  }

  // Aggregate events at same time
  let mut aggregated = BTreeMap::new();
  for &(time, delta) in events {
    *aggregated.entry(OrderedFloat(time)).or_insert(0) += delta;
  }

  // Build breakpoints and values
  // Value before first event is 0 (no lineages yet)
  let mut current_count = 0_i32;
  let (breakpoints, values): (Vec<_>, Vec<_>) = aggregated
    .into_iter()
    .map(|(time, delta)| {
      current_count += delta;
      (time.into_inner(), current_count as f64)
    })
    .unzip();

  let breakpoints = Array1::from_vec(breakpoints);
  let values = std::iter::once(0.0)
    .chain(values)
    .collect_vec();
  let values = Array1::from_vec(values);

  Ok(PiecewiseConstantFn::new(breakpoints, values))
}

#[cfg(test)]
mod tests {
  use ndarray::array;

  use super::*;
  use crate::pretty_assert_abs_diff_eq;

  #[test]
  fn test_piecewise_constant_eval() {
    // Breakpoints: [1.0, 5.0, 10.0]
    // Values: [0.0, 1.0, 2.0, 3.0]
    // t < 1.0 -> 0.0
    // 1.0 <= t < 5.0 -> 1.0
    // 5.0 <= t < 10.0 -> 2.0
    // t >= 10.0 -> 3.0
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_abs_diff_eq!(pc.eval(-1.0), 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(0.5), 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(1.0), 1.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(3.0), 1.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(5.0), 2.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(7.5), 2.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(10.0), 3.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(pc.eval(100.0), 3.0, epsilon = 1e-9);
  }

  #[test]
  fn test_piecewise_constant_eval_many() {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0], array![0.0, 1.0, 2.0]);
    let ts = array![0.0, 1.0, 3.0, 5.0, 10.0];
    let result = pc.eval_many(&ts);
    pretty_assert_abs_diff_eq!(result[0], 0.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(result[1], 1.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(result[2], 1.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(result[3], 2.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(result[4], 2.0, epsilon = 1e-9);
  }

  #[test]
  fn test_lineage_count_simple_tree() -> Result<(), Report> {
    // Simple tree: 2 tips at 0, root at 10.
    // Events:
    // 0.0: Tip 1 (+1)
    // 0.0: Tip 2 (+1)
    // 10.0: Root (-1) (merger of 2 -> 1)
    let events = vec![(0.0, 1), (0.0, 1), (10.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Check values at specific points
    // Before any events (t < 0.0): 0
    pretty_assert_abs_diff_eq!(lineage_counts.eval(-1.0), 0.0, epsilon = 1e-9);
    // After t=0 events (+2): 2
    pretty_assert_abs_diff_eq!(lineage_counts.eval(5.0), 2.0, epsilon = 1e-9);
    // After t=10 event (-1): 1
    pretty_assert_abs_diff_eq!(lineage_counts.eval(15.0), 1.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_single_event() -> Result<(), Report> {
    let events = vec![(5.0, 2)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Before event (t < 5.0): 0
    pretty_assert_abs_diff_eq!(lineage_counts.eval(4.0), 0.0, epsilon = 1e-9);
    // After event (t >= 5.0): 2
    pretty_assert_abs_diff_eq!(lineage_counts.eval(6.0), 2.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_empty_events() {
    let events: Vec<(f64, i32)> = vec![];
    drop(compute_lineage_count_distribution(&events).unwrap_err());
  }

  #[test]
  fn test_lineage_count_aggregation() -> Result<(), Report> {
    let events = vec![(5.0, 1), (5.0, 1), (5.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // Net delta +1 at 5.0
    // Before 5.0: 0
    pretty_assert_abs_diff_eq!(lineage_counts.eval(4.5), 0.0, epsilon = 1e-9);
    // After 5.0: 1
    pretty_assert_abs_diff_eq!(lineage_counts.eval(5.5), 1.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_decreasing_then_increasing() -> Result<(), Report> {
    let events = vec![(0.0, 3), (5.0, -1), (10.0, 1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // 0.0: +3 -> 3
    // 5.0: -1 -> 2
    // 10.0: +1 -> 3
    // Before first event: 0
    pretty_assert_abs_diff_eq!(lineage_counts.eval(-1.0), 0.0, epsilon = 1e-9);
    // After 0.0, before 5.0: 3
    pretty_assert_abs_diff_eq!(lineage_counts.eval(2.5), 3.0, epsilon = 1e-9);
    // After 5.0, before 10.0: 2
    pretty_assert_abs_diff_eq!(lineage_counts.eval(7.5), 2.0, epsilon = 1e-9);
    // After 10.0: 3
    pretty_assert_abs_diff_eq!(lineage_counts.eval(15.0), 3.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_breakpoints() -> Result<(), Report> {
    let events = vec![(10.0, 1), (20.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let breakpoints = lineage_counts.breakpoints();
    assert_eq!(breakpoints.len(), 2);
    pretty_assert_abs_diff_eq!(breakpoints[0], 10.0, epsilon = 1e-9);
    pretty_assert_abs_diff_eq!(breakpoints[1], 20.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_negative_deltas() -> Result<(), Report> {
    let events = vec![(0.0, 5), (5.0, -1), (10.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    // 0.0: +5 -> 5
    // 5.0: -1 -> 4
    // 10.0: -1 -> 3
    // Before first event: 0
    pretty_assert_abs_diff_eq!(lineage_counts.eval(-1.0), 0.0, epsilon = 1e-9);
    // After 0.0, before 5.0: 5
    pretty_assert_abs_diff_eq!(lineage_counts.eval(2.5), 5.0, epsilon = 1e-9);
    // After 5.0, before 10.0: 4
    pretty_assert_abs_diff_eq!(lineage_counts.eval(7.5), 4.0, epsilon = 1e-9);
    // After 10.0: 3
    pretty_assert_abs_diff_eq!(lineage_counts.eval(15.0), 3.0, epsilon = 1e-9);

    Ok(())
  }
}
