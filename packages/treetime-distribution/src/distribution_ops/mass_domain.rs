use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
use crate::distribution_ops::time_bounds::distribution_support_n_points;
use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, s};
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridEdge, GridFn, MaxGridPoints, Side, SoftTailLaw};
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
  let y = f.y();
  let lower = lower_cut(&t, y, f.left_extrap(), fractions.quantile(eps));
  let upper = upper_cut(&t, y, f.right_extrap(), fractions.survival_inverse(eps));
  if lower == MassCut::Kept && upper == MassCut::Kept {
    return Ok(dist);
  }
  let last = t.len() - 1;
  let first = lower.index().unwrap_or(0);
  let end = upper.index().unwrap_or(last).max(first + 1).min(last);
  let first = first.min(end - 1);
  let (first, end) = hard_cell_fit_hull(f.grid_fn(), first, end)?;
  let lower = edge_cut(lower, first, 0, f.left_extrap());
  let upper = edge_cut(upper, end, last, f.right_extrap());
  let left = cut_side(&profile, f.left_extrap(), Side::Left, lower, &t, y)?;
  let right = cut_side(&profile, f.right_extrap(), Side::Right, upper, &t, y)?;
  let trimmed = DistributionFunction::from_start_dx_values(t[first], f.dx(), y.slice(s![first..=end]).to_owned())?
    .with_amplitude(f.amplitude())
    .with_left_extrap(left)?
    .with_right_extrap(right)?;
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

fn lower_cut(t: &Array1<f64>, y: &Array1<f64>, extrap: BoundaryBehavior, quantile: f64) -> MassCut {
  let at_or_below = partition_point(t, |&x| x <= quantile);
  let cut = at_or_below.checked_sub(1).and_then(|index| {
    if extrap.is_soft() {
      (0..=index).rev().find(|&i| y[i].is_finite())
    } else {
      Some(index)
    }
  });
  match cut {
    None => MassCut::Kept,
    Some(0) if extrap.is_soft() => MassCut::Kept,
    Some(index) => MassCut::At(index),
  }
}

fn upper_cut(t: &Array1<f64>, y: &Array1<f64>, extrap: BoundaryBehavior, quantile: f64) -> MassCut {
  let last = t.len() - 1;
  let at_or_above = partition_point(t, |&x| x < quantile);
  let cut = if extrap.is_soft() {
    (at_or_above..=last).find(|&i| y[i].is_finite())
  } else {
    (at_or_above <= last).then_some(at_or_above)
  };
  match cut {
    None => MassCut::Kept,
    Some(index) if index == last && extrap.is_soft() => MassCut::Kept,
    Some(index) => MassCut::At(index),
  }
}

fn hard_cell_fit_hull(grid_fn: &GridFn<f64>, first: usize, end: usize) -> Result<(usize, usize), Report> {
  (first..end).try_fold((first, end), |(lo, hi), cell| {
    Ok(match grid_fn.hard_cell_fit_range(cell)? {
      Some(range) => (lo.min(*range.start()), hi.max(*range.end())),
      None => (lo, hi),
    })
  })
}

fn edge_cut(cut: MassCut, index: usize, edge: usize, extrap: BoundaryBehavior) -> MassCut {
  if index == edge && (cut == MassCut::Kept || extrap.is_soft()) {
    MassCut::Kept
  } else {
    MassCut::At(index)
  }
}

fn cut_side(
  profile: &MassProfile,
  extrap: BoundaryBehavior,
  side: Side,
  cut: MassCut,
  t: &Array1<f64>,
  y: &Array1<f64>,
) -> Result<BoundaryBehavior, Report> {
  let MassCut::At(index) = cut else {
    return Ok(extrap);
  };
  let edge = GridEdge {
    t: t[index],
    y: y[index],
  };
  Ok(match extrap {
    BoundaryBehavior::Linear(_) => BoundaryBehavior::Linear(tail_with_mass_beyond(profile, side, edge)?),
    BoundaryBehavior::HardApproach(law) if ulps_eq!(law.t_hard, edge.t, max_ulps = 4) => extrap,
    BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) | BoundaryBehavior::Error => BoundaryBehavior::Hard,
  })
}

fn tail_with_mass_beyond(profile: &MassProfile, side: Side, edge: GridEdge) -> Result<SoftTailLaw, Report> {
  let neg_log_beyond = match side {
    Side::Left => profile.neg_log_mass_between(f64::NEG_INFINITY, edge.t)?,
    Side::Right => profile.neg_log_mass_between(edge.t, f64::INFINITY)?,
  };
  let rate = (neg_log_beyond - (edge.y - profile.y_peak())).exp();
  let slope = match side {
    Side::Left => -rate,
    Side::Right => rate,
  };
  Ok(SoftTailLaw { slope })
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
