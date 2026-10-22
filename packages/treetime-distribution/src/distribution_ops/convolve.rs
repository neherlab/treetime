use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::mass_profile::MassProfile;
use crate::distribution_ops::multiply::{HardDomain, distribution_hard_domain, guarded_empty_result};
use crate::distribution_ops::time_bounds::distribution_support_n_points;
use crate::policy::NegLog;
use eyre::Report;
use itertools::Itertools;
use ndarray::{Array1, s};
use treetime_grid::{BoundaryBehavior, HardApproachLaw, MaxGridPoints, Side};
use treetime_ops::convolve_fft;
use treetime_utils::array::ndarray::{max_or, min_or};
use treetime_utils::cast::ToF64Exact;
use treetime_utils::{make_error, make_internal_error};

const CONV_TRUST_FRACTION: f64 = 1e-13;

const CONV_TAIL_MARGIN: usize = 3;

const CELLS_PER_RISE: f64 = 8.0;

pub(crate) fn distribution_convolution_fine(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  match (a, b) {
    (Distribution::Function(a), Distribution::Function(b)) => convolution_function_function_fine(a, b, max_points),
    _ => distribution_convolution(a, b, max_points),
  }
}

pub(crate) fn distribution_convolution(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  match (a, b) {
    (Distribution::Formula(_), _)
    | (_, Distribution::Formula(_))
    | (Distribution::Range(_), Distribution::Range(_)) => {
      make_error!("Cannot convolve {a} with {b}: operation not implemented")
    },
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      let a_domain = distribution_hard_domain(a);
      let b_domain = distribution_hard_domain(b);
      guarded_empty_result("convolution", a_domain, b_domain)
    },
    (Distribution::Point(a), Distribution::Point(b)) => Ok(convolution_point_point(a, b)),
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      Ok(convolution_point_range(a, b))
    },
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      Ok(Distribution::Function(convolution_point_function(a, b)?))
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      convolution_range_function(a, b, max_points)
    },
    (Distribution::Function(a), Distribution::Function(b)) => convolution_function_function(a, b, max_points),
  }
}

fn convolution_point_point(
  a: &DistributionPoint<f64, NegLog>,
  b: &DistributionPoint<f64, NegLog>,
) -> Distribution<NegLog> {
  Distribution::point(a.t() + b.t(), a.amplitude() + b.amplitude())
}

fn convolution_point_range(
  p: &DistributionPoint<f64, NegLog>,
  r: &DistributionRange<f64, NegLog>,
) -> Distribution<NegLog> {
  Distribution::range((r.start() + p.t(), r.end() + p.t()), p.amplitude() + r.amplitude())
}

fn convolution_range_function(
  r: &DistributionRange<f64, NegLog>,
  f: &DistributionFunction<f64, NegLog>,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  if f.is_empty() || !f.grid_fn().y_min().is_finite() || r.amplitude() == f64::INFINITY {
    return guarded_empty_result("convolution", conv_operand_domain(false), conv_operand_domain(false));
  }
  if f.len() == 1 {
    let point = DistributionPoint::new(f.x_min(), f.values()[0]);
    return Ok(convolution_point_range(&point, r));
  }
  let left = SupportEnd::of(f, Side::Left);
  let right = SupportEnd::of(f, Side::Right);
  let Some(narrowest) = narrowest_support_component(f, left, right) else {
    return guarded_empty_result("convolution", conv_operand_domain(false), conv_operand_domain(false));
  };
  let profile = MassProfile::new(f)?;
  let (lo, hi) = (left.t + r.start(), right.t + r.end());
  let hard_ends = usize::from(!left.soft) + usize::from(!right.soft);
  let has_hard_cells = hard_ends > 0 || f.y().iter().any(|y| y.is_infinite());
  let rise = (r.end() - r.start()).min(narrowest);
  let dx = if has_hard_cells {
    f.dx().min(rise / CELLS_PER_RISE)
  } else {
    f.dx()
  };
  let n_points = distribution_support_n_points((lo, hi), dx, max_points)?.max(2 + hard_ends);
  let dx = (hi - lo) / (n_points - 1).to_f64_exact();
  let y = (0..n_points)
    .map(|i| {
      let t = if i == n_points - 1 {
        hi
      } else {
        lo + dx * i.to_f64_exact()
      };
      let s_lo = left.clamp_lower(t - r.end());
      let s_hi = right.clamp_upper(t - r.start());
      if s_lo < s_hi {
        Ok(profile.y_peak() + profile.neg_log_mass_between(s_lo, s_hi)?)
      } else {
        Ok(f64::INFINITY)
      }
    })
    .collect::<Result<Array1<f64>, Report>>()?;
  if !y.iter().any(|y| y.is_finite()) {
    return make_internal_error!(
      "Convolution of the range [{}, {}] with a function on [{}, {}] left no point with density at spacing {dx}",
      r.start(),
      r.end(),
      f.x_min(),
      f.x_max()
    );
  }
  let g = DistributionFunction::from_range_values((lo, hi), y)?
    .with_left_extrap(left.result_tail())?
    .with_right_extrap(right.result_tail())?;
  Ok(Distribution::Function(g.with_amplitude(r.amplitude() + f.amplitude())))
}

fn narrowest_support_component(
  f: &DistributionFunction<f64, NegLog>,
  left: SupportEnd,
  right: SupportEnd,
) -> Option<f64> {
  let last = f.len() - 1;
  f.y()
    .iter()
    .enumerate()
    .chunk_by(|(_, y)| y.is_finite())
    .into_iter()
    .filter(|(finite, _)| *finite)
    .filter_map(|(_, run)| {
      let indices = run.map(|(i, _)| i).collect_vec();
      let (&first, &end) = (indices.first()?, indices.last()?);
      let start = match (first, left.soft) {
        (0, true) => f64::NEG_INFINITY,
        (0, false) => left.t,
        _ => grid_x(f, first - 1),
      };
      let end = match (end == last, right.soft) {
        (true, true) => f64::INFINITY,
        (true, false) => right.t,
        (false, _) => grid_x(f, end + 1),
      };
      let width = end - start;
      (width > 0.0).then_some(width)
    })
    .reduce(f64::min)
}

fn grid_x(f: &DistributionFunction<f64, NegLog>, i: usize) -> f64 {
  if i == f.len() - 1 {
    f.x_max()
  } else {
    f.x_min() + f.dx() * i.to_f64_exact()
  }
}

#[derive(Clone, Copy)]
struct SupportEnd {
  behavior: BoundaryBehavior,
  t: f64,
  soft: bool,
}

impl SupportEnd {
  fn of(f: &DistributionFunction<f64, NegLog>, side: Side) -> Self {
    let y = f.y();
    let last = y.len() - 1;
    let (behavior, edge_index, outermost_finite) = match side {
      Side::Left => (f.left_extrap(), 0, y.iter().position(|y| y.is_finite())),
      Side::Right => (f.right_extrap(), last, y.iter().rposition(|y| y.is_finite())),
    };
    match (outermost_finite, side) {
      (Some(i), Side::Left) if i != edge_index => Self::new(BoundaryBehavior::Hard, grid_x(f, i - 1)),
      (Some(i), Side::Right) if i != edge_index => Self::new(BoundaryBehavior::Hard, grid_x(f, i + 1)),
      _ => Self::new(behavior, grid_x(f, edge_index)),
    }
  }

  fn new(behavior: BoundaryBehavior, grid_edge: f64) -> Self {
    let t = match behavior {
      BoundaryBehavior::HardApproach(law) => law.t_hard,
      BoundaryBehavior::Error
      | BoundaryBehavior::Hard
      | BoundaryBehavior::HardEdgeCell(_)
      | BoundaryBehavior::Linear(_) => grid_edge,
    };
    Self {
      behavior,
      t,
      soft: behavior.is_soft(),
    }
  }

  fn clamp_lower(self, s: f64) -> f64 {
    if self.soft { s } else { s.max(self.t) }
  }

  fn clamp_upper(self, s: f64) -> f64 {
    if self.soft { s } else { s.min(self.t) }
  }

  fn result_tail(self) -> BoundaryBehavior {
    match self.behavior {
      BoundaryBehavior::HardApproach(_) | BoundaryBehavior::HardEdgeCell(_) => BoundaryBehavior::Hard,
      other => other,
    }
  }
}

fn convolution_point_function(
  p: &DistributionPoint<f64, NegLog>,
  f: &DistributionFunction<f64, NegLog>,
) -> Result<DistributionFunction<f64, NegLog>, Report> {
  let shifted = DistributionFunction::from_start_dx_values(f.x_min() + p.t(), f.dx(), f.y().clone())?
    .with_left_extrap(shift_tail_law(f.left_extrap(), p.t()))?
    .with_right_extrap(shift_tail_law(f.right_extrap(), p.t()))?;
  Ok(shifted.with_amplitude(f.amplitude() + p.amplitude()))
}

fn shift_tail_law(behavior: BoundaryBehavior, shift: f64) -> BoundaryBehavior {
  match behavior {
    BoundaryBehavior::HardApproach(law) => BoundaryBehavior::HardApproach(HardApproachLaw {
      t_hard: law.t_hard + shift,
      ..law
    }),
    BoundaryBehavior::HardEdgeCell(law) => BoundaryBehavior::HardEdgeCell(law.shift_arg(shift)),
    other => other,
  }
}

fn convolution_function_function(
  a: &DistributionFunction<f64, NegLog>,
  b: &DistributionFunction<f64, NegLog>,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let coarse_dx = a.dx().max(b.dx());
  match convolution_function_function_fine(a, b, max_points)? {
    Distribution::Function(conv_distr) => coarsen_convolution(conv_distr, coarse_dx, max_points),
    other => Ok(other),
  }
}

fn convolution_function_function_fine(
  a: &DistributionFunction<f64, NegLog>,
  b: &DistributionFunction<f64, NegLog>,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  if a.is_empty() || b.is_empty() {
    let a_domain = conv_operand_domain(!a.is_empty());
    let b_domain = conv_operand_domain(!b.is_empty());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  let dx_a = a.dx();
  let dx_b = b.dx();
  let dx = dx_a.min(dx_b);
  if !(dx.is_finite() && dx > 0.0) {
    return make_error!("Invalid grid spacing detected during convolution: {dx}");
  }

  let a = a.resample_dx(dx, max_points)?;
  let b = b.resample_dx(dx, max_points)?;
  if a.is_empty() || b.is_empty() {
    let a_domain = conv_operand_domain(!a.is_empty());
    let b_domain = conv_operand_domain(!b.is_empty());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  let amplitude = a.amplitude() + b.amplitude();

  let (pa, peak_a) = to_peak_normalized_plain(a.y());
  let (pb, peak_b) = to_peak_normalized_plain(b.y());
  if !(peak_a.is_finite() && peak_b.is_finite()) {
    let a_domain = conv_operand_domain(peak_a.is_finite());
    let b_domain = conv_operand_domain(peak_b.is_finite());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  let result_range = (a.x_min() + b.x_min(), a.x_max() + b.x_max());
  max_points.point_count((a.len() + b.len() - 1).to_f64_exact(), result_range, dx)?;
  let conv = convolve_fft(dx, &pa, &pb)?;

  let Some(reconstructed) = reconstruct_neg_log_tails(&conv, dx)? else {
    let a_domain = conv_operand_domain(true);
    let b_domain = conv_operand_domain(true);
    return guarded_empty_result("convolution", a_domain, b_domain);
  };

  let x_min = a.x_min() + b.x_min() + reconstructed.start_offset.to_f64_exact() * dx;

  let offset = peak_a + peak_b;
  let y = reconstructed.neg_log.mapv(|v| v + offset);

  if y.len() == 1 {
    return Ok(Distribution::point(x_min, amplitude + y[0]));
  }

  let conv_distr = DistributionFunction::from_start_dx_values(x_min, dx, y)?.with_amplitude(amplitude);
  Ok(Distribution::Function(conv_distr))
}

fn conv_operand_domain(has_mass: bool) -> Option<HardDomain> {
  has_mass.then_some((
    (f64::NEG_INFINITY, f64::INFINITY),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
  ))
}

pub(super) fn coarsen_convolution(
  conv_distr: DistributionFunction<f64, NegLog>,
  coarse_dx: f64,
  max_points: MaxGridPoints,
) -> Result<Distribution<NegLog>, Report> {
  let range = conv_distr.x_max() - conv_distr.x_min();
  let coarse_intervals = (range / coarse_dx).round();
  if coarse_intervals.is_nan() || coarse_intervals < 1.0 {
    return Ok(Distribution::Function(conv_distr));
  }
  Ok(Distribution::Function(conv_distr.resample_dx(coarse_dx, max_points)?))
}

fn to_peak_normalized_plain(neg_log: &Array1<f64>) -> (Array1<f64>, f64) {
  let peak = min_or(neg_log, f64::INFINITY);
  if !peak.is_finite() {
    return (Array1::zeros(neg_log.len()), peak);
  }
  let plain = neg_log.mapv(|nl| (peak - nl).exp());
  (plain, peak)
}

#[expect(clippy::integer_division, reason = "integer division is the intended floor division")]
fn reconstruct_neg_log_tails(conv: &Array1<f64>, dx: f64) -> Result<Option<ReconstructedConv>, Report> {
  let peak = max_or(conv, 0.0);
  if peak <= 0.0 {
    return Ok(None);
  }
  let floor = peak * CONV_TRUST_FRACTION;
  let (Some(first), Some(last)) = (
    conv.iter().position(|&c| c > floor),
    conv.iter().rposition(|&c| c > floor),
  ) else {
    return Ok(None);
  };

  let trusted_y = conv.slice(s![first..=last]).mapv(|c| -c.ln());
  let n = trusted_y.len();
  let margin = CONV_TAIL_MARGIN.min(n / 3);
  if margin < 1 {
    return make_error!("Convolution left too few trusted points to reconstruct tails");
  }

  let left_slope = (trusted_y[margin] - trusted_y[0]) / (margin.to_f64_exact() * dx);
  let right_slope = (trusted_y[n - 1] - trusted_y[n - 1 - margin]) / (margin.to_f64_exact() * dx);

  let left_start = if first > 0 && left_slope < 0.0 { 0 } else { first };
  let right_end = if last + 1 < conv.len() && right_slope > 0.0 {
    conv.len() - 1
  } else {
    last
  };

  let mut neg_log = Array1::<f64>::zeros(right_end - left_start + 1);
  for i in first..=last {
    neg_log[i - left_start] = trusted_y[i - first];
  }
  for i in left_start..first {
    neg_log[i - left_start] = trusted_y[0] + left_slope * ((i.to_f64_exact() - first.to_f64_exact()) * dx);
  }
  for i in (last + 1)..=right_end {
    neg_log[i - left_start] = trusted_y[n - 1] + right_slope * ((i.to_f64_exact() - last.to_f64_exact()) * dx);
  }

  Ok(Some(ReconstructedConv {
    start_offset: left_start,
    neg_log,
  }))
}

struct ReconstructedConv {
  start_offset: usize,
  neg_log: Array1<f64>,
}
