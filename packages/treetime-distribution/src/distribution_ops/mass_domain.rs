use crate::distribution_ops::mass_profile::{MassProfile, TailMass};
use crate::distribution_ops::time_bounds::distribution_support_n_points;
use crate::policy::NegLog;
use crate::{Distribution, DistributionFunction};
use approx::ulps_eq;
use eyre::Report;
use itertools::Either;
use ndarray::{Array1, s};
use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, GridEdge, GridFn, MaxGridPoints, Side, SoftTailLaw};
use treetime_utils::array::ndarray::partition_point;
use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};
use treetime_utils::{make_error, make_internal_error};

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
  let cuts = TrimCuts {
    profile: &profile,
    grid_fn: f.grid_fn(),
    t: &t,
    y,
    left: f.left_extrap(),
    right: f.right_extrap(),
  };
  let lower = cuts.lower(fractions.quantile(eps))?;
  let upper = cuts.upper(fractions.survival_inverse(eps))?;
  if lower == MassCut::Kept && upper == MassCut::Kept {
    return Ok(dist);
  }
  let window = cuts.settled_window(lower, upper)?;
  let left = cuts.cut_side(Side::Left, window.lower)?;
  let right = cuts.cut_side(Side::Right, window.upper)?;
  let trimmed = DistributionFunction::from_start_dx_values(
    t[window.first],
    f.dx(),
    y.slice(s![window.first..=window.end]).to_owned(),
  )?
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

struct TrimCuts<'a> {
  profile: &'a MassProfile,
  grid_fn: &'a GridFn<f64>,
  t: &'a Array1<f64>,
  y: &'a Array1<f64>,
  left: BoundaryBehavior,
  right: BoundaryBehavior,
}

impl TrimCuts<'_> {
  fn lower(&self, quantile: f64) -> Result<MassCut, Report> {
    let at_or_below = partition_point(self.t, |&x| x <= quantile);
    match at_or_below.checked_sub(1) {
      None => Ok(MassCut::Kept),
      Some(index) => self.outward(Side::Left, index),
    }
  }

  fn upper(&self, quantile: f64) -> Result<MassCut, Report> {
    let at_or_above = partition_point(self.t, |&x| x < quantile);
    if at_or_above > self.last() {
      return Ok(MassCut::Kept);
    }
    self.outward(Side::Right, at_or_above)
  }

  fn settled_window(&self, lower: MassCut, upper: MassCut) -> Result<TrimWindow, Report> {
    let mut cuts = (lower, upper);
    for _ in 0..=2 * self.t.len() {
      let window = self.window(cuts.0, cuts.1)?;
      let moved = (
        self.recut(Side::Left, window.lower)?,
        self.recut(Side::Right, window.upper)?,
      );
      if moved == (window.lower, window.upper) {
        return Ok(window);
      }
      cuts = moved;
    }
    make_internal_error!(
      "Mass trim cuts did not settle within the grid of {} points",
      self.t.len()
    )
  }

  fn window(&self, lower: MassCut, upper: MassCut) -> Result<TrimWindow, Report> {
    let last = self.last();
    let first = lower.index().unwrap_or(0);
    let end = upper.index().unwrap_or(last).max(first + 1).min(last);
    let first = first.min(end - 1);
    let (first, end) = hard_cell_fit_hull(self.grid_fn, first, end)?;
    Ok(TrimWindow {
      first,
      end,
      lower: edge_cut(lower, first, 0, self.left),
      upper: edge_cut(upper, end, last, self.right),
    })
  }

  fn recut(&self, side: Side, cut: MassCut) -> Result<MassCut, Report> {
    match cut {
      MassCut::At(index) => self.outward(side, index),
      MassCut::Kept => Ok(MassCut::Kept),
    }
  }

  fn outward(&self, side: Side, index: usize) -> Result<MassCut, Report> {
    if !self.extrap(side).is_soft() {
      return Ok(MassCut::At(index));
    }
    let candidates = match side {
      Side::Left => Either::Left((1..=index).rev()),
      Side::Right => Either::Right(index..self.last()),
    };
    let found = candidates
      .map(|i| self.soft_tail_at(side, i).map(|law| law.map(|_| i)))
      .find_map(Result::transpose)
      .transpose()?;
    Ok(found.map_or(MassCut::Kept, MassCut::At))
  }

  fn cut_side(&self, side: Side, cut: MassCut) -> Result<BoundaryBehavior, Report> {
    let extrap = self.extrap(side);
    let MassCut::At(index) = cut else {
      return Ok(extrap);
    };
    Ok(match extrap {
      BoundaryBehavior::Linear(_) => match self.soft_tail_at(side, index)? {
        Some(law) => BoundaryBehavior::Linear(law),
        None => return make_internal_error!("Mass trim cut at grid point {index} has no representable tail"),
      },
      BoundaryBehavior::HardApproach(law) if ulps_eq!(law.t_hard, self.t[index], max_ulps = 4) => extrap,
      BoundaryBehavior::Hard | BoundaryBehavior::HardApproach(_) | BoundaryBehavior::Error => BoundaryBehavior::Hard,
    })
  }

  fn soft_tail_at(&self, side: Side, index: usize) -> Result<Option<SoftTailLaw>, Report> {
    let edge = GridEdge {
      t: self.t[index],
      y: self.y[index],
    };
    if !edge.y.is_finite() {
      return Ok(None);
    }
    tail_with_mass_beyond(self.profile, side, edge)
  }

  fn extrap(&self, side: Side) -> BoundaryBehavior {
    match side {
      Side::Left => self.left,
      Side::Right => self.right,
    }
  }

  fn last(&self) -> usize {
    self.t.len() - 1
  }
}

struct TrimWindow {
  first: usize,
  end: usize,
  lower: MassCut,
  upper: MassCut,
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

fn tail_with_mass_beyond(profile: &MassProfile, side: Side, edge: GridEdge) -> Result<Option<SoftTailLaw>, Report> {
  let neg_log_beyond = match side {
    Side::Left => profile.neg_log_mass_between(f64::NEG_INFINITY, edge.t)?,
    Side::Right => profile.neg_log_mass_between(edge.t, f64::INFINITY)?,
  };
  let rate = (neg_log_beyond - (edge.y - profile.y_peak())).exp();
  if !rate.is_normal() {
    return Ok(None);
  }
  let slope = match side {
    Side::Left => -rate,
    Side::Right => rate,
  };
  Ok(Some(SoftTailLaw { slope }))
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
