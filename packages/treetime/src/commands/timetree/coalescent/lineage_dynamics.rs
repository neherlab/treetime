use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use eyre::Report;
use ndarray::{Array1, s};
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use treetime_utils::make_error;

/// Computes k(t) distribution from tree events.
///
/// k(t) is the number of concurrent lineages at time t.
/// The function is piecewise constant, stepping at each merger event.
/// Events must be sorted by increasing time (past to present).
pub fn compute_lineage_count_distribution(events: &[(f64, i32)]) -> Result<Distribution, Report> {
  if events.is_empty() {
    return make_error!("Cannot build lineage count interpolator from empty events");
  }

  let mut aggregated_events = BTreeMap::new();
  for &(time, delta) in events {
    *aggregated_events.entry(OrderedFloat(time)).or_insert(0) += delta;
  }

  let mut cumulative_events: Vec<(f64, i32)> = Vec::with_capacity(aggregated_events.len());
  let mut current_count = 0_i32;
  // Iterate from present (t=0) to past (t>0)
  for (time, delta) in aggregated_events {
    cumulative_events.push((time.into_inner(), current_count));
    current_count += delta;
  }

  let data_t_min = cumulative_events.first().unwrap().0;
  let data_t_max = cumulative_events.last().unwrap().0;
  let data_range = data_t_max - data_t_min;
  let margin = f64::max(data_range * 0.1, 1.0);

  let t_min = data_t_min - margin;
  let t_max = data_t_max + margin;
  let n_points = 10000;
  let dx = (t_max - t_min) / (n_points - 1) as f64;

  let mut y_vals = Array1::zeros(n_points);
  let mut prev_grid_idx = 0;

  for (event_time, count) in cumulative_events {
    let next_grid_idx = (((event_time - t_min) / dx).ceil() as usize).min(n_points);

    if prev_grid_idx < next_grid_idx {
      y_vals.slice_mut(s![prev_grid_idx..next_grid_idx]).fill(count as f64);
    }

    prev_grid_idx = next_grid_idx;
  }

  // Fill the remaining part (past the last event) with the final count
  if prev_grid_idx < n_points {
    y_vals.slice_mut(s![prev_grid_idx..n_points]).fill(current_count as f64);
  }

  let func = DistributionFunction::from_range_values((t_min, t_max), y_vals)?;
  Ok(Distribution::Function(func))
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::pretty_assert_abs_diff_eq;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_lineage_count_simple_tree() -> Result<(), Report> {
    // Simple tree: 2 tips at 0, root at 10.
    // Events:
    // 0.0: Tip 1 (+1)
    // 0.0: Tip 2 (+1)
    // 10.0: Root (-1) (merger of 2 -> 1)
    let events = vec![(0.0, 1), (0.0, 1), (10.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    assert_eq!(t_grid.len(), 10000);

    let idx_past = t_grid.iter().position(|&t| t < -0.5).unwrap();
    let idx_mid = t_grid.iter().position(|&t| t > 2.0 && t < 4.0).unwrap();
    let idx_high = t_grid.iter().position(|&t| t > 10.5).unwrap();

    // Before tips (future): 0
    pretty_assert_abs_diff_eq!(y_grid[idx_past], 0.0, epsilon = 1e-9);
    // Between tips and root: 2
    pretty_assert_abs_diff_eq!(y_grid[idx_mid], 2.0, epsilon = 1e-9);
    // After root (past): 1
    pretty_assert_abs_diff_eq!(y_grid[idx_high], 1.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_single_event() -> Result<(), Report> {
    let events = vec![(5.0, 2)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    assert_eq!(t_grid.len(), 10000);

    let idx_before = t_grid.iter().position(|&t| (4.0..5.0).contains(&t)).unwrap();
    let idx_after = t_grid.iter().position(|&t| t >= 6.0).unwrap();

    // Before event (present side): 0
    pretty_assert_abs_diff_eq!(y_grid[idx_before], 0.0, epsilon = 1e-9);
    // After event (past side): 2
    pretty_assert_abs_diff_eq!(y_grid[idx_after], 2.0, epsilon = 1e-9);

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

    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    let idx_before = t_grid.iter().position(|&t| t > 4.0 && t < 4.9).unwrap();
    let idx_after = t_grid.iter().position(|&t| t > 5.1 && t < 6.0).unwrap();

    // Net delta +1 at 5.0
    // Before 5.0: 0
    pretty_assert_abs_diff_eq!(y_grid[idx_before], 0.0, epsilon = 1e-9);
    // After 5.0: 1
    pretty_assert_abs_diff_eq!(y_grid[idx_after], 1.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_decreasing_then_increasing() -> Result<(), Report> {
    let events = vec![(0.0, 3), (5.0, -1), (10.0, 1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    let idx_mid = t_grid.iter().position(|&t| t > 6.0 && t < 9.0).unwrap();

    // 0.0: +3 -> 3
    // 5.0: -1 -> 2
    // 10.0: +1 -> 3
    // Mid (6..9): 2
    pretty_assert_abs_diff_eq!(y_grid[idx_mid], 2.0, epsilon = 1e-9);

    Ok(())
  }

  #[test]
  fn test_lineage_count_margin_bounds() -> Result<(), Report> {
    let events = vec![(10.0, 1), (20.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let t_grid = lineage_counts.t();

    let expected_margin = 1.0;

    assert!(t_grid[0] <= 10.0 - expected_margin + 0.1);
    assert!(t_grid[t_grid.len() - 1] >= 20.0 + expected_margin - 0.1);

    Ok(())
  }

  #[test]
  fn test_lineage_count_negative_deltas() -> Result<(), Report> {
    let events = vec![(0.0, 5), (5.0, -1), (10.0, -1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    let idx_mid = t_grid.iter().position(|&t| t > 6.0 && t < 9.0).unwrap();

    // 0.0: +5 -> 5
    // 5.0: -1 -> 4
    // 10.0: -1 -> 3
    // Mid (6..9): 4
    pretty_assert_abs_diff_eq!(y_grid[idx_mid], 4.0, epsilon = 1e-9);

    Ok(())
  }
}
