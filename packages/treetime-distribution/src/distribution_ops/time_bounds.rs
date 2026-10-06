use eyre::Report;
use treetime_grid::MaxGridPoints;
use treetime_utils::make_error;

const GRID_COUNT_TOL: f64 = 1e-9;

pub(super) fn distribution_support_n_points(
  (start, end): (f64, f64),
  dx: f64,
  max_points: MaxGridPoints,
) -> Result<usize, Report> {
  let ratio = (end - start) / dx;
  let nearest = ratio.round();
  let intervals = if (ratio - nearest).abs() <= GRID_COUNT_TOL * ratio.abs().max(1.0) {
    nearest
  } else {
    ratio.ceil()
  };
  if intervals.is_nan() || intervals < 0.0 {
    return make_error!("Cannot discretize distribution support [{start}, {end}] with spacing {dx}");
  }
  let n_points = max_points.point_count(intervals + 1.0, (start, end), dx)?;
  Ok(n_points.max(2))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum SupportIntersection {
  Disjoint,
  Point(f64),
  Interval((f64, f64)),
}
