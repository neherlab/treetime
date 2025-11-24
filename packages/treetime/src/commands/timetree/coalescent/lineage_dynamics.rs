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
  for (time, delta) in aggregated_events.into_iter().rev() {
    current_count += delta;
    cumulative_events.push((time.into_inner(), current_count));
  }
  cumulative_events.reverse();

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

  let func = DistributionFunction::from_range_values((t_min, t_max), y_vals)?;
  Ok(Distribution::Function(func))
}

#[cfg(test)]
mod tests {
  use super::*;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_lineage_count_simple_tree() -> Result<(), Report> {
    let events = vec![(10.0, 1), (10.0, 1), (5.0, -1), (0.0, 1)];

    let lineage_counts = compute_lineage_count_distribution(&events)?;

    let t_grid = lineage_counts.t();
    let y_grid = lineage_counts.y();

    assert!(t_grid.len() == 10000, "Expected 10000 grid points");

    let idx_future = t_grid.iter().position(|&t| t > 10.5).unwrap();
    let idx_high = t_grid.iter().position(|&t| t > 7.0 && t < 9.0).unwrap();
    let idx_mid = t_grid.iter().position(|&t| t > 2.0 && t < 4.0).unwrap();
    let idx_past = t_grid.iter().position(|&t| t < -0.5).unwrap();

    assert_abs_diff_eq!(y_grid[idx_future], 0.0, epsilon = 0.1);
    assert_abs_diff_eq!(y_grid[idx_high], 2.0, epsilon = 0.1);
    assert_abs_diff_eq!(y_grid[idx_mid], 1.0, epsilon = 0.1);
    assert_abs_diff_eq!(y_grid[idx_past], 2.0, epsilon = 0.1);

    Ok(())
  }
}
