use crate::Distribution;
use crate::policy::YAxisPolicy;
use eyre::Report;
use num::ToPrimitive;
use treetime_utils::make_error;

pub(super) const MAX_GRID_POINTS: usize = 1_000_000;

/// Compute union of time bounds from two distributions.
///
/// Return (t_min, t_max) tuple representing the smallest interval that encompasses
/// both distributions' time domains.
pub fn distribution_time_bounds_union<Y: YAxisPolicy>(
  dist_a: &Distribution<Y>,
  dist_b: &Distribution<Y>,
) -> (f64, f64) {
  let (t_min_a, t_max_a) = dist_a.time_bounds();
  let (t_min_b, t_max_b) = dist_b.time_bounds();
  let t_min = f64::min(t_min_a, t_min_b);
  let t_max = f64::max(t_max_a, t_max_b);
  (t_min, t_max)
}

/// Compute intersection of time bounds from two distributions.
///
/// Return Some((t_min, t_max)) representing the overlapping time interval
/// or None if intervals don't overlap.
pub fn distribution_time_bounds_intersection<Y: YAxisPolicy>(
  dist_a: &Distribution<Y>,
  dist_b: &Distribution<Y>,
) -> Option<(f64, f64)> {
  match distribution_support_intersection(dist_a.time_bounds(), dist_b.time_bounds()) {
    SupportIntersection::Disjoint => None,
    SupportIntersection::Point(t) => Some((t, t)),
    SupportIntersection::Interval(bounds) => Some(bounds),
  }
}

/// Whether inner distribution's time bounds are fully contained within outer distribution's time bounds.
pub fn distribution_time_bounds_contains<Y: YAxisPolicy>(outer: &Distribution<Y>, inner: &Distribution<Y>) -> bool {
  let (t_min_outer, t_max_outer) = outer.time_bounds();
  let (t_min_inner, t_max_inner) = inner.time_bounds();
  t_min_outer <= t_min_inner && t_max_inner <= t_max_outer
}

/// Check if two distributions' time bounds overlap.
pub fn distribution_time_bounds_overlaps<Y: YAxisPolicy>(dist_a: &Distribution<Y>, dist_b: &Distribution<Y>) -> bool {
  distribution_time_bounds_intersection(dist_a, dist_b).is_some()
}

#[allow(clippy::float_cmp)] // Point support requires exact contact; a tolerance would enlarge the intersection.
pub(super) fn distribution_support_intersection(a: (f64, f64), b: (f64, f64)) -> SupportIntersection {
  let start = a.0.max(b.0);
  let end = a.1.min(b.1);

  if start < end {
    SupportIntersection::Interval((start, end))
  } else if start == end {
    SupportIntersection::Point(start)
  } else {
    SupportIntersection::Disjoint
  }
}

pub(super) fn distribution_support_n_points((start, end): (f64, f64), dx: f64) -> Result<usize, Report> {
  let intervals = ((end - start) / dx).round();
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
