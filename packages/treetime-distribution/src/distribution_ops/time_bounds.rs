use crate::Distribution;
use crate::policy::YAxisPolicy;
use eyre::Report;
use num::ToPrimitive;
use treetime_utils::make_error;

pub(super) const MAX_GRID_POINTS: usize = 1_000_000;

/// Relative tolerance for treating an interval-count ratio as an exact integer when sizing a grid.
/// Far above the few-ulp error of a single division, far below a half interval, so it absorbs float
/// noise at exact multiples without ever misclassifying a genuine fractional ratio.
const GRID_COUNT_TOL: f64 = 1e-9;

/// Compute union of time bounds from two distributions.
///
/// `Empty` is the union identity: it contributes nothing. Two empty distributions
/// yield `None`.
pub fn distribution_time_bounds_union<Y: YAxisPolicy>(
  dist_a: &Distribution<Y>,
  dist_b: &Distribution<Y>,
) -> Option<(f64, f64)> {
  match (dist_a.time_bounds(), dist_b.time_bounds()) {
    (Some((t_min_a, t_max_a)), Some((t_min_b, t_max_b))) => {
      Some((f64::min(t_min_a, t_min_b), f64::max(t_max_a, t_max_b)))
    },
    (Some(bounds), None) | (None, Some(bounds)) => Some(bounds),
    (None, None) => None,
  }
}

/// Compute intersection of time bounds from two distributions.
///
/// Return Some((t_min, t_max)) representing the overlapping time interval
/// or None if intervals don't overlap. `Empty` is absorbing: intersecting with
/// it yields `None`.
pub fn distribution_time_bounds_intersection<Y: YAxisPolicy>(
  dist_a: &Distribution<Y>,
  dist_b: &Distribution<Y>,
) -> Option<(f64, f64)> {
  match distribution_support_intersection(dist_a.time_bounds()?, dist_b.time_bounds()?) {
    SupportIntersection::Disjoint => None,
    SupportIntersection::Point(t) => Some((t, t)),
    SupportIntersection::Interval(bounds) => Some(bounds),
  }
}

/// Whether inner distribution's time bounds are fully contained within outer distribution's
/// time bounds. The empty support is contained in any distribution; no non-empty support is
/// contained in an empty one.
pub fn distribution_time_bounds_contains<Y: YAxisPolicy>(outer: &Distribution<Y>, inner: &Distribution<Y>) -> bool {
  let Some((t_min_inner, t_max_inner)) = inner.time_bounds() else {
    return true;
  };
  let Some((t_min_outer, t_max_outer)) = outer.time_bounds() else {
    return false;
  };
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
  // Round the interval count up (resolution floor): the realized spacing (end - start) / intervals is
  // then never coarser than the finest operand's dx, while `Array1::linspace` still lands both
  // analytical endpoints on grid nodes. See kb/decisions/distribution-tails-and-arithmetic.md.
  //
  // A ratio within a small relative tolerance of an integer is taken as that integer, so float error
  // at an exact multiple (e.g. 1.5 / 0.5 evaluating to 3.0000000000000004) does not add a spurious
  // point. The tolerance sits far above the few-ulp division error and far below a half interval.
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
