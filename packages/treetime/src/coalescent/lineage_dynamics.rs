use crate::coalescent::time_coordinate::CalendarTime;
use eyre::Report;
use itertools::Itertools;
use ndarray::Array1;
use ordered_float::OrderedFloat;
use std::collections::BTreeMap;
use std::iter::once;
use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
use treetime_utils::make_error;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn compute_lineage_count_distribution(
  events: &[(CalendarTime, i32)],
  terminal_lineage_count: i32,
) -> Result<PiecewiseConstantFn, Report> {
  if events.is_empty() {
    return make_error!("Cannot build lineage count from empty events");
  }
  if terminal_lineage_count < 0 {
    return make_error!("Terminal lineage count must be non-negative, got {terminal_lineage_count}");
  }

  let mut aggregated = BTreeMap::new();
  for &(time, delta) in events {
    *aggregated.entry(OrderedFloat(time.value())).or_insert(0) += delta;
  }

  let mut current_count = 1_i32;
  let (breakpoints, values): (Vec<_>, Vec<_>) = aggregated
    .into_iter()
    .map(|(time, delta)| {
      current_count -= delta;
      (time.into_inner(), current_count as f64)
    })
    .unzip();

  if current_count != terminal_lineage_count {
    return make_error!(
      "Lineage count must end at {terminal_lineage_count} after the latest retained sample, got {current_count}"
    );
  }

  let breakpoints = Array1::from_vec(breakpoints);
  let values = once(1.0).chain(values).collect_vec();
  let values = Array1::from_vec(values);

  Ok(PiecewiseConstantFn::new(breakpoints, values))
}
