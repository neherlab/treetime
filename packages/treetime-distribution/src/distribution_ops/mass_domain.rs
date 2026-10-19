use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
use crate::distribution_ops::time_bounds::distribution_support_n_points;
use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::make_error;

pub(crate) fn log_mass(dist: &Distribution<NegLog>) -> Result<f64, Report> {
  match dist {
    Distribution::Empty => Ok(f64::NEG_INFINITY),
    Distribution::Point(p) => Ok(-p.amplitude()),
    Distribution::Range(r) => Ok(-r.amplitude() + (r.end() - r.start()).ln()),
    Distribution::Function(f) => {
      let profile = MassProfile::new(f)?;
      let relative = profile.neg_log_mass_between(f64::NEG_INFINITY, f64::INFINITY)?;
      Ok(-f.amplitude() - profile.y_peak() - relative)
    },
    Distribution::Formula(_) => log_mass(&dist.repeak()?),
  }
}

pub fn rewindow_to_mass(
  dist: &Distribution<NegLog>,
  eps: f64,
  grid_points: usize,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let Distribution::Function(f) = dist else {
    return dist.repeak();
  };
  let Some(sized) = peak_normalized_if_mass_sizable(f) else {
    return dist.repeak();
  };
  let (lo, hi) = sized.bounded_domain(eps)?;
  resample_to_mass_window(&sized.normalized, lo, hi, grid_points, max_points)
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    error_dropped_by_pattern,
    reason = "a mass profile that cannot be computed means the distribution is not mass-sizable, and callers fall back to plain normalization"
  )
)]
pub(crate) fn peak_normalized_if_mass_sizable(f: &DistributionFunction<f64, NegLog>) -> Option<MassSized> {
  let Ok(Some(normalized)) = f.repeak() else {
    return None;
  };
  let Ok(profile) = MassProfile::new(&normalized) else {
    return None;
  };
  is_mass_sizable(&profile).then_some(MassSized { normalized, profile })
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    error_dropped_by_pattern,
    reason = "a mass profile that cannot be computed means the distribution is not mass-sizable, and callers fall back to plain normalization"
  )
)]
pub(crate) fn mass_sizable_profile(f: &DistributionFunction<f64, NegLog>) -> Option<MassProfile> {
  let Ok(profile) = MassProfile::new(f) else {
    return None;
  };
  is_mass_sizable(&profile).then_some(profile)
}

pub(crate) struct MassSized {
  pub(crate) normalized: DistributionFunction<f64, NegLog>,
  profile: MassProfile,
}

impl MassSized {
  pub(crate) fn bounded_domain(&self, eps: f64) -> Result<(f64, f64), Report> {
    mass_bounded_domain(&self.profile, eps)
  }
}

pub(crate) fn mass_bounded_domain(profile: &MassProfile, eps: f64) -> Result<(f64, f64), Report> {
  let Some(fractions) = profile.fractions(TailMass::Required) else {
    return make_error!(
      "Mass-bounded domain requires a finite positive total mass; the distribution is not mass-sizable"
    );
  };
  let lo = match profile.extrap(Side::Left) {
    BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) => profile.x_min(),
    BoundaryBehavior::Linear(_) => fractions.quantile(eps),
    BoundaryBehavior::Error => return make_error!("Mass-bounded domain: left side has a non-mass-sizable tail"),
  };
  let hi = match profile.extrap(Side::Right) {
    BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) => profile.x_max(),
    BoundaryBehavior::Linear(_) => fractions.survival_inverse(eps),
    BoundaryBehavior::Error => return make_error!("Mass-bounded domain: right side has a non-mass-sizable tail"),
  };
  if hi <= lo {
    return make_error!("Mass-bounded domain collapsed to an empty interval [{lo}, {hi}]");
  }
  Ok((lo, hi))
}

pub(crate) fn resample_to_mass_window(
  normalized: &DistributionFunction<f64, NegLog>,
  lo: f64,
  hi: f64,
  grid_points: usize,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  if hi <= lo {
    return make_error!("Mass window collapsed to an empty interval [{lo}, {hi}]");
  }
  let mass_width = hi - lo;
  let floor_dx = normalized.dx();
  let target_dx = mass_width / grid_points.saturating_sub(1).max(1).to_f64_exact();
  #[expect(
    clippy::float_cmp,
    reason = "a window bounded by the hard edges of the input grid is that grid, and only identical edges give its point count"
  )]
  let window_is_input_grid = lo == normalized.x_min() && hi == normalized.x_max();
  let n_points = if target_dx <= floor_dx {
    max_points.point_count(grid_points.max(2).to_f64_exact(), (lo, hi), target_dx)?
  } else if window_is_input_grid {
    max_points.point_count(normalized.len().to_f64_exact(), (lo, hi), floor_dx)?
  } else {
    distribution_support_n_points((lo, hi), floor_dx, max_points)?
  };
  let resampled = normalized.resample_range_n_points((lo, hi), n_points)?;
  let resampled = refit_soft_tails(resampled)?;
  Ok(Distribution::Function(resampled))
}

pub(crate) fn refit_soft_tails(
  f: DistributionFunction<f64, NegLog>,
) -> Result<DistributionFunction<f64, NegLog>, Report> {
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

fn is_mass_sizable(profile: &MassProfile) -> bool {
  profile.fractions(TailMass::Required).is_some()
}
