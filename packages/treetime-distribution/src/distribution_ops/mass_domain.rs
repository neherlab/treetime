use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
use crate::distribution_ops::time_bounds::distribution_support_n_points;
use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use eyre::Report;
use ndarray::{Array1, s};
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridFn, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::array::ndarray::partition_point;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
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
  resample_to_mass_window(&sized.normalized, lo, hi, grid_points, max_points, SpacingLimit::Input)
}

pub fn trim_to_mass(dist: Distribution<NegLog>, eps: f64) -> Result<Distribution<NegLog>, Report> {
  let Distribution::Function(f) = &dist else {
    return Ok(dist);
  };
  let Some(profile) = mass_sizable_profile(f) else {
    return Ok(dist);
  };
  let Some(fractions) = profile.fractions(TailMass::Required) else {
    return Ok(dist);
  };
  let t = f.t();
  let lower = lower_cut(&t, f.left_extrap(), fractions.quantile(eps));
  let upper = upper_cut(&t, f.right_extrap(), fractions.survival_inverse(eps));
  if lower == MassCut::Kept && upper == MassCut::Kept {
    return Ok(dist);
  }
  let last = t.len() - 1;
  let first = lower.index().unwrap_or(0);
  let end = upper.index().unwrap_or(last).max(first + 1).min(last);
  let first = first.min(end - 1);
  let trimmed = DistributionFunction::from_start_dx_values(t[first], f.dx(), f.y().slice(s![first..=end]).to_owned())?
    .with_amplitude(f.amplitude())
    .with_left_extrap(cut_side(f.left_extrap(), lower))?
    .with_right_extrap(cut_side(f.right_extrap(), upper))?;
  let trimmed = match lower {
    MassCut::At(_) if trimmed.left_extrap().is_soft() => {
      let law = SoftTailLaw::fit(trimmed.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS)?;
      trimmed.with_left_extrap(BoundaryBehavior::Linear(law))?
    },
    MassCut::At(_) | MassCut::Kept => trimmed,
  };
  let trimmed = match upper {
    MassCut::At(_) if trimmed.right_extrap().is_soft() => {
      let law = SoftTailLaw::fit(trimmed.grid_fn(), Side::Right, DEFAULT_TAIL_FIT_POINTS)?;
      trimmed.with_right_extrap(BoundaryBehavior::Linear(law))?
    },
    MassCut::At(_) | MassCut::Kept => trimmed,
  };
  Ok(Distribution::Function(trimmed))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MassCut {
  Kept,
  At(usize),
}

impl MassCut {
  fn index(self) -> Option<usize> {
    match self {
      MassCut::Kept => None,
      MassCut::At(index) => Some(index),
    }
  }
}

fn lower_cut(t: &Array1<f64>, extrap: BoundaryBehavior, quantile: f64) -> MassCut {
  let at_or_below = partition_point(t, |&x| x <= quantile);
  match at_or_below.checked_sub(1) {
    None => MassCut::Kept,
    Some(0) if extrap.is_soft() => MassCut::Kept,
    Some(index) => MassCut::At(index),
  }
}

fn upper_cut(t: &Array1<f64>, extrap: BoundaryBehavior, quantile: f64) -> MassCut {
  let last = t.len() - 1;
  match partition_point(t, |&x| x < quantile) {
    index if index > last => MassCut::Kept,
    index if index == last && extrap.is_soft() => MassCut::Kept,
    index => MassCut::At(index),
  }
}

fn cut_side(extrap: BoundaryBehavior, cut: MassCut) -> BoundaryBehavior {
  match cut {
    MassCut::At(_) if !extrap.is_soft() => BoundaryBehavior::Hard,
    MassCut::At(_) | MassCut::Kept => extrap,
  }
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
  spacing: SpacingLimit,
) -> Result<Distribution<NegLog>, Report> {
  if hi <= lo {
    return make_error!("Mass window collapsed to an empty interval [{lo}, {hi}]");
  }
  let mass_width = hi - lo;
  let max_dx = match spacing {
    SpacingLimit::Input => normalized.dx(),
    SpacingLimit::AtMost(max_dx) => max_dx,
  };
  let target_dx = mass_width / grid_points.saturating_sub(1).max(1).to_f64_exact();
  #[expect(
    clippy::float_cmp,
    reason = "a window bounded by the hard edges of the input grid is that grid, and only identical edges give its point count"
  )]
  let window_is_input_grid =
    matches!(spacing, SpacingLimit::Input) && lo == normalized.x_min() && hi == normalized.x_max();
  let n_points = if target_dx <= max_dx {
    max_points.point_count(grid_points.max(2).to_f64_exact(), (lo, hi), target_dx)?
  } else if window_is_input_grid {
    max_points.point_count(normalized.len().to_f64_exact(), (lo, hi), max_dx)?
  } else {
    distribution_support_n_points((lo, hi), max_dx, max_points)?
  };
  let resampled = normalized.resample_range_n_points((lo, hi), n_points)?;
  let resampled = match spacing {
    SpacingLimit::Input => refit_soft_tails(resampled)?,
    SpacingLimit::AtMost(_) => soft_tails_at_window_edges(normalized, resampled, (lo, hi))?,
  };
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

fn soft_tails_at_window_edges(
  input: &DistributionFunction<f64, NegLog>,
  resampled: DistributionFunction<f64, NegLog>,
  (lo, hi): (f64, f64),
) -> Result<DistributionFunction<f64, NegLog>, Report> {
  let resampled = if matches!(resampled.left_extrap(), BoundaryBehavior::Linear(_)) {
    let law = window_edge_tail_law(input, &resampled, Side::Left, lo)?;
    resampled.with_left_extrap(BoundaryBehavior::Linear(law))?
  } else {
    resampled
  };
  let resampled = if matches!(resampled.right_extrap(), BoundaryBehavior::Linear(_)) {
    let law = window_edge_tail_law(input, &resampled, Side::Right, hi)?;
    resampled.with_right_extrap(BoundaryBehavior::Linear(law))?
  } else {
    resampled
  };
  Ok(resampled)
}

fn window_edge_tail_law(
  input: &DistributionFunction<f64, NegLog>,
  resampled: &DistributionFunction<f64, NegLog>,
  side: Side,
  edge: f64,
) -> Result<SoftTailLaw, Report> {
  const NODE_ROUNDING: f64 = 1e-9;
  let last = input.len() - 1;
  let offset = (edge - input.x_min()) / input.dx();
  let (first, end) = match side {
    Side::Left => {
      let first = f64_to_usize_saturating((offset - NODE_ROUNDING).ceil()).min(last);
      (first, (first + DEFAULT_TAIL_FIT_POINTS - 1).min(last))
    },
    Side::Right => {
      let end = f64_to_usize_saturating((offset + NODE_ROUNDING).floor()).min(last);
      ((end + 1).saturating_sub(DEFAULT_TAIL_FIT_POINTS), end)
    },
  };
  let finite = input
    .y()
    .slice(s![first..=end])
    .iter()
    .filter(|y| y.is_finite())
    .count();
  if finite < 2 {
    return SoftTailLaw::fit(resampled.grid_fn(), side, DEFAULT_TAIL_FIT_POINTS);
  }
  let window_edge = GridFn::from_start_dx_values(
    input.x_min() + first.to_f64_exact() * input.dx(),
    input.dx(),
    input.y().slice(s![first..=end]).to_owned(),
  )?;
  SoftTailLaw::fit(&window_edge, side, DEFAULT_TAIL_FIT_POINTS)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum SpacingLimit {
  Input,
  AtMost(f64),
}

fn is_mass_sizable(profile: &MassProfile) -> bool {
  profile.fractions(TailMass::Required).is_some()
}
