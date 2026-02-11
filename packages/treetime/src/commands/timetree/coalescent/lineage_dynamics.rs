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
  let values = std::iter::once(0.0).chain(values).collect_vec();
  let values = Array1::from_vec(values);

  Ok(PiecewiseConstantFn::new(breakpoints, values))
}
