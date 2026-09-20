use crate::policy::{NegLog, YAxisPolicy};
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::Array1;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridEdge, Side, SoftTailLaw};
use treetime_utils::make_error;

pub fn total_mass(f: &DistributionFunction<f64, NegLog>) -> Result<f64, Report> {
  Ok(mass_profile(f)?.z)
}

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

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
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

fn soft_edge(profile: &MassProfile, slope: f64, tail_mass: f64, target: f64, side: Side) -> f64 {
  if tail_mass >= target {
    let shift = (tail_mass / target).ln() / slope;
    match side {
      Side::Left => profile.x_min - shift,
      Side::Right => profile.x_max + shift,
    }
  } else {
    trim_into_grid(profile, target - tail_mass, side)
  }
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn trim_into_grid(profile: &MassProfile, grid_target: f64, side: Side) -> f64 {
  let plain = &profile.plain;
  let dx = profile.dx;
  let n_points = plain.len();

  let mut acc = 0.0;
  for step in 0..n_points - 1 {
    let j = match side {
      Side::Left => step,
      Side::Right => n_points - 2 - step,
    };
    let dp = plain[j + 1] - plain[j];
    let cell_mass = 0.5 * (plain[j] + plain[j + 1]) * dx;
    if acc + cell_mass >= grid_target {
      let remaining = (grid_target - acc) / dx;
      let q = match side {
        Side::Left => remaining,
        Side::Right => plain[j] + 0.5 * dp - remaining,
      };
      let h = solve_cell_fraction(plain[j], dp, q);
      return profile.x_min + (j as f64 + h) * dx;
    }
    acc += cell_mass;
  }
  match side {
    Side::Left => profile.x_max,
    Side::Right => profile.x_min,
  }
}

fn solve_cell_fraction(p0: f64, dp: f64, q: f64) -> f64 {
  let disc = (p0 * p0 + 2.0 * dp * q).max(0.0);
  let denom = p0 + disc.sqrt();
  let h = if denom > 0.0 { 2.0 * q / denom } else { 0.0 };
  h.clamp(0.0, 1.0)
}

struct MassProfile {
  plain: Array1<f64>,
  dx: f64,
  x_min: f64,
  x_max: f64,
  left: BoundaryBehavior,
  right: BoundaryBehavior,
  left_mass: f64,
  right_mass: f64,
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

fn tail_mass(extrap: BoundaryBehavior, w_edge: f64, t_edge: f64, side: Side) -> Result<f64, Report> {
  match extrap {
    BoundaryBehavior::Linear(law) => Ok(law.mass(w_edge)),
    BoundaryBehavior::HardApproach(law) => Ok(law.mass(GridEdge { t: t_edge, y: w_edge })),
    BoundaryBehavior::Hard => Ok(0.0),
    BoundaryBehavior::Error => make_error!("Mass domain: {side:?} side has an undeclared (Error) tail"),
  }
}
