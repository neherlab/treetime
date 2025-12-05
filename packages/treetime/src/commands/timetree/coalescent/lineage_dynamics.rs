use eyre::Report;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_utils::make_error;

/// Piecewise constant function represented by breakpoints and values.
///
/// For n breakpoints, there are n+1 value regions:
/// - values[0] for t < breakpoints[0]
/// - values[i] for breakpoints[i-1] <= t < breakpoints[i]  (1 <= i < n)
/// - values[n] for t >= breakpoints[n-1]
///
/// Breakpoints must be sorted in ascending order.
#[derive(Debug, Clone)]
pub struct PiecewiseConstant {
  breakpoints: Vec<f64>,
  values: Vec<f64>,
}

impl PiecewiseConstant {
  /// Create from breakpoints and values.
  ///
  /// # Arguments
  /// - `breakpoints`: sorted ascending, n elements
  /// - `values`: n+1 elements
  pub fn new(breakpoints: Vec<f64>, values: Vec<f64>) -> Self {
    debug_assert!(breakpoints.len() + 1 == values.len());
    debug_assert!(breakpoints
      .windows(2)
      .all(|w| matches!(w, [a, b] if a < b)));
    Self { breakpoints, values }
  }

  /// Get breakpoint times (where function changes value).
  pub fn breakpoints(&self) -> &[f64] {
    &self.breakpoints
  }

  /// Evaluate at a single point.
  pub fn eval(&self, t: f64) -> f64 {
    let idx = self.breakpoints.partition_point(|&bp| bp <= t);
    self.values[idx]
  }

  /// Evaluate at multiple points.
  pub fn eval_many(&self, ts: &Array1<f64>) -> Array1<f64> {
    ts.mapv(|t| self.eval(t))
  }
}

/// Computes k(t) distribution from tree events.
///
/// k(t) is the number of concurrent lineages at time t.
/// The function is piecewise constant, stepping at each merger event.
/// Events must be sorted by increasing time (past to present).
pub fn compute_lineage_count_distribution(events: &[(f64, i32)]) -> Result<PiecewiseConstant, Report> {
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
  let mut breakpoints = Vec::with_capacity(aggregated.len());
  let mut values = Vec::with_capacity(aggregated.len() + 1);
  values.push(0.0); // value for t < first_breakpoint

  let mut current_count = 0_i32;
  for (time, delta) in aggregated {
    breakpoints.push(time.into_inner());
    current_count += delta;
    values.push(current_count as f64);
  }

  Ok(PiecewiseConstant::new(breakpoints, values))
}

#[cfg(test)]
mod tests {
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
    let pc = PiecewiseConstant::new(vec![1.0, 5.0, 10.0], vec![0.0, 1.0, 2.0, 3.0]);

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
    let pc = PiecewiseConstant::new(vec![1.0, 5.0], vec![0.0, 1.0, 2.0]);
    let ts = Array1::from_vec(vec![0.0, 1.0, 3.0, 5.0, 10.0]);
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
