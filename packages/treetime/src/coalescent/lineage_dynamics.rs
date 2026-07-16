use crate::coalescent::time_coordinate::CalendarTime;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::iter::once;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::make_error;

/// Computes k(t) distribution from tree events in calendar-year coordinates.
///
/// k(t) is the number of concurrent lineages at time t.
/// The function is piecewise constant, stepping at each merger event.
/// Events must be sorted by increasing time (past to present). Event deltas are
/// expressed in the time-before-present direction, so they are subtracted while
/// traversing calendar time toward the present.
pub fn compute_lineage_count_distribution(events: &[(CalendarTime, i32)]) -> Result<PiecewiseConstantFn, Report> {
  if events.is_empty() {
    return make_error!("Cannot build lineage count from empty events");
  }

  // Aggregate events at same time
  let mut aggregated = BTreeMap::new();
  for &(time, delta) in events {
    *aggregated.entry(OrderedFloat(time.value())).or_insert(0) += delta;
  }

  // Older than the root, the sampled tree has one ancestral lineage. Moving
  // toward the present reverses the event deltas collected for TBP traversal.
  let mut current_count = 1_i32;
  let (breakpoints, values): (Vec<_>, Vec<_>) = aggregated
    .into_iter()
    .map(|(time, delta)| {
      current_count -= delta;
      (time.into_inner(), current_count as f64)
    })
    .unzip();

  if current_count != 0 {
    return make_error!("Lineage count must end at zero after the latest sample, got {current_count}");
  }

  let breakpoints = Array1::from_vec(breakpoints);
  let values = once(1.0).chain(values).collect_vec();
  let values = Array1::from_vec(values);

  Ok(PiecewiseConstantFn::new(breakpoints, values))
}
