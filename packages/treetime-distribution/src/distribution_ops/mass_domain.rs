use crate::policy::{NegLog, YAxisPolicy};
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridEdge, Side, SoftTailLaw};
use treetime_utils::make_error;

/// Total probability mass of a neg-log grid distribution, tails included.
///
/// The mass is measured peak-relative: the peak ordinate is subtracted before exponentiating, so the
/// mode has plain probability one and the grid trapezoid and both tail masses share one reference. An
/// absolute (non-peak-relative) edge ordinate would mis-scale the total by `exp(y_peak)` because the
/// closed-form tail masses ([`SoftTailLaw::mass`], [`HardApproachLaw::mass`](treetime_grid::HardApproachLaw::mass))
/// read the same stored edge ordinate.
///
/// `Z = trapezoid(exp(-(y - y_peak))) + left_tail_mass + right_tail_mass`, where each tail mass is the
/// closed-form half-line (soft) or sub-grid-gap (hard-approach) integral, and a `Hard` edge
/// contributes nothing (probability is exactly zero beyond it). A degenerate soft tail with zero decay
/// (a flat, non-integrable likelihood) contributes infinite mass, so the total is `+inf`; that is a
/// value, not an error, and lets the re-window fall back gracefully. Returns an error only for a side
/// whose tail is genuinely non-mass-sizable and never appears on a stored inference distribution
/// (an undeclared `Error` tail).
pub fn total_mass(f: &DistributionFunction<f64, NegLog>) -> Result<f64, Report> {
  Ok(mass_profile(f)?.z)
}

/// The mass-bounded grid domain `(lo, hi)` for a neg-log grid distribution.
///
/// Each side is sized independently by its declared tail (design D2, soft-side-only trim):
///
/// - soft (`Linear`): the edge is the point where the mass outside it (its tail plus the grid mass
///   integrated inward) equals `eps * Z`. Solved in closed form within the tail when the tail alone
///   already carries more than `eps * Z`, and by trapezoid-CDF inversion within the grid otherwise.
/// - `Hard`: the grid edge is the hard bound itself, so it stays put; there is no mass beyond it. A
///   finite mode sitting exactly on the bound is carried here, on the exact grid endpoint.
/// - `HardApproach` (divergent): the grid edge is strictly inside `t_hard` with a sub-grid gap, so no
///   `+inf` ordinate is ever placed on the grid; the edge stays at the grid bound.
///
/// Probability is exactly zero beyond a hard bound, so `eps` never trims a hard side. Returns an error
/// for an `Error`-tailed side.
pub fn mass_bounded_domain(f: &DistributionFunction<f64, NegLog>, eps: f64) -> Result<(f64, f64), Report> {
  let profile = mass_profile(f)?;
  if !(profile.z.is_finite() && profile.z > 0.0) {
    return make_error!(
      "Mass-bounded domain requires a finite positive total mass, got {}; the distribution is not mass-sizable",
      profile.z
    );
  }
  let target = eps * profile.z;
  let lo = lower_edge(&profile, target)?;
  let hi = upper_edge(&profile, target)?;
  if hi <= lo {
    return make_error!("Mass-bounded domain collapsed to an empty interval [{lo}, {hi}]");
  }
  Ok((lo, hi))
}

/// Peak-normalize a neg-log distribution and resample it over its mass-bounded domain.
///
/// This is the one re-window rule the timetree passes apply after every composite step. It first
/// subtracts the peak ordinate (subsuming the shift-only `normalize` it replaces), computes the mass
/// domain (D2), and resamples at the resolution floor
/// `dx = min(mass_width / (grid_points - 1), input_dx)`, so the result holds **at least**
/// `grid_points` points and is never coarser than the input. [`GridFn::resample`](treetime_grid::GridFn)
/// carries both edge-relative tail laws across the domain change, so a hard side and the raw soft law
/// survive automatically; the soft tails are then re-fit on the moved outermost points as a deliberate
/// accuracy step at the new edge.
///
/// A non-`Function` distribution (`Point`, `Range`, `Empty`, `Formula`) has no grid to re-window, so it
/// returns the shift-only `normalize`, exactly what this call replaces at those sites.
pub fn rewindow_to_mass(
  dist: &Distribution<NegLog>,
  eps: f64,
  grid_points: usize,
) -> Result<Distribution<NegLog>, Report> {
  let Distribution::Function(f) = dist else {
    return Ok(dist.normalize());
  };
  let Some(normalized) = peak_normalized_if_mass_sizable(f) else {
    return Ok(dist.normalize());
  };
  let (lo, hi) = mass_bounded_domain(&normalized, eps)?;
  resample_to_mass_window(&normalized, lo, hi, grid_points)
}

/// Peak-normalize a neg-log grid distribution, or `None` when it is not mass-sizable.
///
/// Subtracts the peak ordinate (so the mode has plain probability one) when the distribution has a
/// finite mode and a finite positive total mass. A distribution with no finite mode, a degenerate soft
/// tail with zero decay (infinite mass), or an `Error` tail carrying no probability law -- such as the
/// hard box a leaf date range contributes to a forward product -- is not mass-sizable; the caller
/// falls back to the shift-only `normalize`, since mass-sizing is an accuracy step, not a correctness
/// requirement.
pub fn peak_normalized_if_mass_sizable(
  f: &DistributionFunction<f64, NegLog>,
) -> Option<DistributionFunction<f64, NegLog>> {
  let y_peak = f.grid_fn().y_min();
  if !y_peak.is_finite() {
    return None;
  }
  let normalized = f.shift_y(-y_peak);
  total_mass(&normalized)
    .is_ok_and(|z| z.is_finite() && z > 0.0)
    .then_some(normalized)
}

/// Resample a peak-normalized neg-log distribution onto the mass window `[lo, hi]`.
///
/// The resolution floor is `dx = min(mass_width / (grid_points - 1), input_dx)`, so the result holds
/// **at least** `grid_points` points and is never coarser than the input. The resample runs over the
/// exact `[lo, hi]` endpoints rather than a `dx`-stepped grid whose final point can round a fraction of
/// a cell past the bound: a hard edge must land exactly on its bound, so a downstream division by a
/// hard-bounded message never samples it beyond support (which would read `+inf` and collapse the
/// message to empty). The `ceil` keeps the resulting spacing `mass_width / (n - 1) <= dx`, so the
/// floor holds. The soft tails are then re-fit on the moved outermost points as a deliberate accuracy
/// step at the new edge.
pub fn resample_to_mass_window(
  normalized: &DistributionFunction<f64, NegLog>,
  lo: f64,
  hi: f64,
  grid_points: usize,
) -> Result<Distribution<NegLog>, Report> {
  if hi <= lo {
    return make_error!("Mass window collapsed to an empty interval [{lo}, {hi}]");
  }
  let mass_width = hi - lo;
  let floor_dx = normalized.dx();
  let target_dx = mass_width / (grid_points.saturating_sub(1).max(1) as f64);
  let dx = target_dx.min(floor_dx);
  let n_points = ((mass_width / dx).ceil() as usize + 1).max(grid_points).max(2);
  let resampled = normalized.resample_range_n_points((lo, hi), n_points)?;
  let resampled = refit_soft_tails(resampled)?;
  Ok(Distribution::Function(resampled))
}

/// Refit the soft (`Linear`) tail on each moved edge from the outermost grid points.
///
/// Only soft tails are refit: their slope is recovered by least squares from the grid, so moving the
/// edge (or changing the ordinates by a non-constant amount) is a genuine accuracy improvement. A
/// `HardApproach` law needs the producer-supplied boundary ordinate that the grid does not carry, so
/// it is left unchanged (edge-relative and already valid); `Hard` and `Error` carry no fittable
/// slope and are left unchanged too.
///
/// Shared with [`distribution_multiply_by_fn`](crate::distribution_multiply_by_fn), which adds
/// a per-point cost to the ordinates and re-fits the soft tail whose slope the cost changed.
pub fn refit_soft_tails(f: DistributionFunction<f64, NegLog>) -> Result<DistributionFunction<f64, NegLog>, Report> {
  let f = if matches!(f.left_extrap(), BoundaryBehavior::Linear(_)) {
    let law = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?;
    f.with_left_extrap(BoundaryBehavior::Linear(law))?
  } else {
    f
  };
  let f = if matches!(f.right_extrap(), BoundaryBehavior::Linear(_)) {
    let law = SoftTailLaw::fit(f.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
    f.with_right_extrap(BoundaryBehavior::Linear(law))?
  } else {
    f
  };
  Ok(f)
}

/// Lower (left) mass-domain edge for the profile's declared left tail (design D2).
fn lower_edge(profile: &MassProfile, target: f64) -> Result<f64, Report> {
  match profile.left {
    BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) => Ok(profile.x_min),
    BoundaryBehavior::Linear(law) => Ok(soft_edge(
      profile,
      law.slope.abs(),
      profile.left_mass,
      target,
      Side::Left,
    )),
    BoundaryBehavior::Error => {
      make_error!("Mass-bounded domain: left side has a non-mass-sizable tail")
    },
  }
}

/// Upper (right) mass-domain edge for the profile's declared right tail (design D2).
fn upper_edge(profile: &MassProfile, target: f64) -> Result<f64, Report> {
  match profile.right {
    BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) => Ok(profile.x_max),
    BoundaryBehavior::Linear(law) => Ok(soft_edge(
      profile,
      law.slope.abs(),
      profile.right_mass,
      target,
      Side::Right,
    )),
    BoundaryBehavior::Error => {
      make_error!("Mass-bounded domain: right side has a non-mass-sizable tail")
    },
  }
}

/// Soft-side edge where the outer mass equals `target = eps * Z`.
///
/// When the tail alone already carries at least `target`, the edge extends into the tail so the
/// remaining outer mass shrinks to `target` (closed form on the exponential tail). Otherwise the edge
/// trims inward so the tail plus the trimmed grid mass reaches `target` (trapezoid-CDF inversion).
fn soft_edge(profile: &MassProfile, slope: f64, tail_mass: f64, target: f64, side: Side) -> f64 {
  if tail_mass >= target {
    // Extend outward: exp(-slope * shift) = target / tail_mass, shift >= 0 since tail_mass >= target.
    let shift = (tail_mass / target).ln() / slope;
    match side {
      Side::Left => profile.x_min - shift,
      Side::Right => profile.x_max + shift,
    }
  } else {
    trim_into_grid(profile, target - tail_mass, side)
  }
}

/// Trim a soft edge inward until the grid mass beyond it equals `grid_target`, by trapezoid-CDF
/// inversion. Accumulates whole-cell masses from the edge inward, then solves the partial cell.
fn trim_into_grid(profile: &MassProfile, grid_target: f64, side: Side) -> f64 {
  let plain = &profile.plain;
  let dx = profile.dx;
  let n_points = plain.len();

  let mut acc = 0.0;
  for step in 0..n_points - 1 {
    // Walk cells inward from the trimmed edge: left-to-right on the left side, right-to-left on the right.
    let j = match side {
      Side::Left => step,
      Side::Right => n_points - 2 - step,
    };
    let dp = plain[j + 1] - plain[j];
    let cell_mass = 0.5 * (plain[j] + plain[j + 1]) * dx;
    if acc + cell_mass >= grid_target {
      let remaining = (grid_target - acc) / dx;
      // Fraction of the cell measured from the outward edge inward.
      let q = match side {
        Side::Left => remaining,
        Side::Right => plain[j] + 0.5 * dp - remaining,
      };
      let h = solve_cell_fraction(plain[j], dp, q);
      return profile.x_min + (j as f64 + h) * dx;
    }
    acc += cell_mass;
  }
  // Unreachable for a valid profile: grid_target < grid_mass by construction. Fall back to the far edge.
  match side {
    Side::Left => profile.x_max,
    Side::Right => profile.x_min,
  }
}

/// Solve `0.5 * dp * h^2 + p0 * h = q` for the fraction `h` in `[0, 1]` within one linear cell.
///
/// Uses the cancellation-free root `h = 2q / (p0 + sqrt(p0^2 + 2 dp q))`, which stays exact as
/// `dp -> 0` (the uniform-cell limit `h = q / p0`) with no division by `dp`.
fn solve_cell_fraction(p0: f64, dp: f64, q: f64) -> f64 {
  let disc = (p0 * p0 + 2.0 * dp * q).max(0.0);
  let denom = p0 + disc.sqrt();
  let h = if denom > 0.0 { 2.0 * q / denom } else { 0.0 };
  h.clamp(0.0, 1.0)
}

/// Peak-relative mass decomposition of a neg-log grid distribution.
struct MassProfile {
  /// Plain-space, peak-relative ordinates `exp(-(y_i - y_peak))` at each grid point.
  plain: Array1<f64>,
  dx: f64,
  x_min: f64,
  x_max: f64,
  left: BoundaryBehavior,
  right: BoundaryBehavior,
  /// Mass in the left tail beyond `x_min`.
  left_mass: f64,
  /// Mass in the right tail beyond `x_max`.
  right_mass: f64,
  /// Total mass `Z` (grid trapezoid plus both tails).
  z: f64,
}

fn mass_profile(f: &DistributionFunction<f64, NegLog>) -> Result<MassProfile, Report> {
  let ys = f.y();
  let n_points = ys.len();
  if n_points < 2 {
    return make_error!("Mass domain needs at least two grid points, got {n_points}");
  }

  let y_peak = f.grid_fn().y_min();
  if !y_peak.is_finite() {
    return make_error!("Mass domain requires a finite peak ordinate, got {y_peak}");
  }

  let plain = ys.mapv(|yi| NegLog::to_plain(yi - y_peak));
  let dx = f.dx();
  let x_min = f.x_min();
  let x_max = f.x_max();
  let grid_mass = dx * (plain.sum() - 0.5 * (plain[0] + plain[n_points - 1]));

  let left = f.left_extrap();
  let right = f.right_extrap();
  let left_mass = tail_mass(left, ys[0] - y_peak, x_min, Side::Left)?;
  let right_mass = tail_mass(right, ys[n_points - 1] - y_peak, x_max, Side::Right)?;

  // `total` may be `+inf` when a soft tail has zero decay (a flat, non-integrable distribution). That
  // is a legitimate value here; `mass_bounded_domain` rejects it, and the re-window falls back.
  let total = grid_mass + left_mass + right_mass;

  Ok(MassProfile {
    plain,
    dx,
    x_min,
    x_max,
    left,
    right,
    left_mass,
    right_mass,
    z: total,
  })
}

/// Closed-form tail mass beyond one grid edge, peak-relative.
///
/// `w_edge` is the peak-relative neg-log edge ordinate (`y_edge - y_peak`) and `t_edge` its coordinate.
/// A `Hard` edge terminates the domain, so its mass is zero. An `Error` edge is not mass-sizable and
/// is rejected.
fn tail_mass(extrap: BoundaryBehavior, w_edge: f64, t_edge: f64, side: Side) -> Result<f64, Report> {
  match extrap {
    BoundaryBehavior::Linear(law) => Ok(law.mass(w_edge)),
    BoundaryBehavior::HardApproach(law) => Ok(law.mass(GridEdge { t: t_edge, y: w_edge })),
    BoundaryBehavior::Hard => Ok(0.0),
    BoundaryBehavior::Error => make_error!("Mass domain: {side:?} side has an undeclared (Error) tail"),
  }
}
