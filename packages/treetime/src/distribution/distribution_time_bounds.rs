use crate::distribution::distribution::Distribution;
use crate::distribution::y_axis_policy::YAxisPolicy;

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
  let (t_min_a, t_max_a) = dist_a.time_bounds();
  let (t_min_b, t_max_b) = dist_b.time_bounds();
  let t_min = f64::max(t_min_a, t_min_b);
  let t_max = f64::min(t_max_a, t_max_b);
  (t_min <= t_max).then_some((t_min, t_max))
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
