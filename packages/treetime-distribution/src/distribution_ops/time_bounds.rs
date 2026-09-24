use eyre::Report;
use num::ToPrimitive;
use treetime_utils::make_error;

const MAX_GRID_POINTS: usize = 1_000_000;

const GRID_COUNT_TOL: f64 = 1e-9;

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(super) fn distribution_support_n_points((start, end): (f64, f64), dx: f64) -> Result<usize, Report> {
  let ratio = (end - start) / dx;
  let nearest = ratio.round();
  let intervals = if (ratio - nearest).abs() <= GRID_COUNT_TOL * ratio.abs().max(1.0) {
    nearest
  } else {
    ratio.ceil()
  };
  if !intervals.is_finite() || intervals < 0.0 {
    return make_error!("Cannot discretize distribution support [{start}, {end}] with spacing {dx}");
  }
  if intervals >= (MAX_GRID_POINTS - 1) as f64 {
    return Ok(MAX_GRID_POINTS);
  }
  let Some(n_points) = intervals.to_usize().and_then(|intervals| intervals.checked_add(1)) else {
    return make_error!("Distribution support [{start}, {end}] with spacing {dx} exceeds the grid size limit");
  };
  Ok(n_points.clamp(2, MAX_GRID_POINTS))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum SupportIntersection {
  Disjoint,
  Point(f64),
  Interval((f64, f64)),
}
