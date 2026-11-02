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
use ndarray::{Array1, Zip};
use std::borrow::Cow;
use std::collections::VecDeque;
use treetime_grid::{BoundaryBehavior, Grid, GridFn, MaxGridPoints, Side};
use treetime_ops::{convolve_direct_at, convolve_fft};
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
      Ok(Distribution::Function(convolution_point_function(a, b)))
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
  let grid = Grid::from_range_n_points(lo, hi, n_points)?;
  let y = grid
    .iter()
    .map(|t| {
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
      "Convolution of the range [{}, {}] with a function on [{}, {}] left no point with density at spacing {}",
      r.start(),
      r.end(),
      f.x_min(),
      f.x_max(),
      grid.dx()
    );
  }
  let g = DistributionFunction::from_grid_fn(GridFn::from_grid_values(grid, y)?)
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
        _ => f.x_at(first - 1),
      };
      let end = match (end == last, right.soft) {
        (true, true) => f64::INFINITY,
        (true, false) => right.t,
        (false, _) => f.x_at(end + 1),
      };
      let width = end - start;
      (width > 0.0).then_some(width)
    })
    .reduce(f64::min)
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
      (Some(i), Side::Left) if i != edge_index => Self::new(BoundaryBehavior::Hard, f.x_at(i - 1)),
      (Some(i), Side::Right) if i != edge_index => Self::new(BoundaryBehavior::Hard, f.x_at(i + 1)),
      _ => Self::new(behavior, f.x_at(edge_index)),
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
) -> DistributionFunction<f64, NegLog> {
  f.shift_time(p.t()).with_amplitude(f.amplitude() + p.amplitude())
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

  let a = on_common_spacing(a, dx_b, max_points)?;
  let b = on_common_spacing(b, dx_a, max_points)?;
  if a.is_empty() || b.is_empty() {
    let a_domain = conv_operand_domain(!a.is_empty());
    let b_domain = conv_operand_domain(!b.is_empty());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  let amplitude = a.amplitude() + b.amplitude();

  let operands = ConvOperands {
    a: PeakNormalizedOperand::new(a.y()),
    b: PeakNormalizedOperand::new(b.y()),
    hard_left: !a.left_extrap().is_soft() && !b.left_extrap().is_soft(),
    hard_right: !a.right_extrap().is_soft() && !b.right_extrap().is_soft(),
  };
  if !(operands.a.peak.is_finite() && operands.b.peak.is_finite()) {
    let a_domain = conv_operand_domain(operands.a.peak.is_finite());
    let b_domain = conv_operand_domain(operands.b.peak.is_finite());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  let result_range = (a.x_min() + b.x_min(), a.x_max() + b.x_max());
  max_points.point_count((a.len() + b.len() - 1).to_f64_exact(), result_range, dx)?;
  let conv = convolve_fft(dx, &operands.a.plain, &operands.b.plain)?;

  let Some(reconstructed) = reconstruct_neg_log_tails(&conv, dx, &operands)? else {
    let a_domain = conv_operand_domain(true);
    let b_domain = conv_operand_domain(true);
    return guarded_empty_result("convolution", a_domain, b_domain);
  };

  let full_grid = Grid::from_ends_dx(result_range.0, result_range.1, dx, conv.len())?;
  let first = reconstructed.start_offset;
  let last = first + reconstructed.neg_log.len() - 1;

  let offset = operands.a.peak + operands.b.peak;
  let y = reconstructed.neg_log.mapv(|v| v + offset);

  if y.len() == 1 {
    return Ok(Distribution::point(full_grid.x_at(first), amplitude + y[0]));
  }

  let grid_fn = GridFn::from_grid_values(full_grid.slice(first..=last)?, y)?;
  let conv_distr = DistributionFunction::from_grid_fn(grid_fn).with_amplitude(amplitude);
  Ok(Distribution::Function(conv_distr))
}

pub(crate) fn on_common_spacing(
  operand: &DistributionFunction<f64, NegLog>,
  other_dx: f64,
  max_points: MaxGridPoints,
) -> Result<Cow<'_, DistributionFunction<f64, NegLog>>, Report> {
  if operand.dx() <= other_dx {
    Ok(Cow::Borrowed(operand))
  } else {
    operand.resample_dx(other_dx, max_points).map(Cow::Owned)
  }
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

#[expect(clippy::integer_division, reason = "integer division is the intended floor division")]
fn reconstruct_neg_log_tails(
  conv: &Array1<f64>,
  dx: f64,
  operands: &ConvOperands,
) -> Result<Option<ReconstructedConv>, Report> {
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

  let trusted_y =
    neg_log_with_quadrature_in_valleys(conv, first, last, floor, dx, &operands.a.plain, &operands.b.plain)?;
  let n = trusted_y.len();
  let margin = CONV_TAIL_MARGIN.min(n / 3);
  if margin < 1 {
    return make_error!("Convolution left too few trusted points to reconstruct tails");
  }

  let trusted_slopes = EndSlopes::of(&trusted_y, margin, dx);
  let left_start = if first > 0 && trusted_slopes.left < 0.0 {
    0
  } else {
    first
  };
  let right_end = if last + 1 < conv.len() && trusted_slopes.right > 0.0 {
    conv.len() - 1
  } else {
    last
  };

  let resolved = neg_log_at_ends_by_tilting(trusted_y, first, left_start, right_end, margin, dx, operands)?;
  let resolved_end = resolved.start + resolved.neg_log.len() - 1;
  let slopes = EndSlopes::of(&resolved.neg_log, margin, dx);
  let left_edge = resolved.neg_log[0];
  let right_edge = resolved.neg_log[resolved.neg_log.len() - 1];
  let start = if operands.hard_left { resolved.start } else { left_start };
  let end = if operands.hard_right { resolved_end } else { right_end };

  let neg_log = (start..=end)
    .map(|i| {
      if i < resolved.start {
        left_edge + slopes.left * ((i.to_f64_exact() - resolved.start.to_f64_exact()) * dx)
      } else if i > resolved_end {
        right_edge + slopes.right * ((i.to_f64_exact() - resolved_end.to_f64_exact()) * dx)
      } else {
        resolved.neg_log[i - resolved.start]
      }
    })
    .collect();

  Ok(Some(ReconstructedConv {
    start_offset: start,
    neg_log,
  }))
}

fn neg_log_with_quadrature_in_valleys(
  conv: &Array1<f64>,
  first: usize,
  last: usize,
  floor: f64,
  dx: f64,
  pa: &Array1<f64>,
  pb: &Array1<f64>,
) -> Result<Array1<f64>, Report> {
  (first..=last)
    .map(|i| {
      let c = if conv[i] > floor {
        conv[i]
      } else {
        convolve_direct_at(dx, pa, pb, i)?.max(f64::MIN_POSITIVE)
      };
      Ok(-c.ln())
    })
    .collect()
}

fn neg_log_at_ends_by_tilting(
  trusted_y: Array1<f64>,
  first: usize,
  left_start: usize,
  right_end: usize,
  margin: usize,
  dx: f64,
  operands: &ConvOperands,
) -> Result<ResolvedNegLog, Report> {
  let mut values: VecDeque<f64> = trusted_y.into_iter().collect();
  let mut start = first;

  while start + values.len() <= right_end {
    let slope = (values[values.len() - 1] - values[values.len() - 1 - margin]) / (margin.to_f64_exact() * dx);
    if slope <= 0.0 {
      break;
    }
    let tilted = TiltedConvolution::new(operands, dx, slope)?;
    let next = start + values.len();
    let added = (next..=right_end).map_while(|i| tilted.neg_log_at(i)).collect_vec();
    if added.is_empty() {
      break;
    }
    values.extend(added);
  }

  while start > left_start {
    let slope = (values[margin] - values[0]) / (margin.to_f64_exact() * dx);
    if slope >= 0.0 {
      break;
    }
    let tilted = TiltedConvolution::new(operands, dx, slope)?;
    let added = (left_start..start)
      .rev()
      .map_while(|i| tilted.neg_log_at(i))
      .collect_vec();
    if added.is_empty() {
      break;
    }
    start -= added.len();
    added.into_iter().for_each(|value| values.push_front(value));
  }

  Ok(ResolvedNegLog {
    start,
    neg_log: values.into_iter().collect(),
  })
}

struct ConvOperands {
  a: PeakNormalizedOperand,
  b: PeakNormalizedOperand,
  hard_left: bool,
  hard_right: bool,
}

struct PeakNormalizedOperand {
  peak: f64,
  relative_neg_log: Array1<f64>,
  plain: Array1<f64>,
}

impl PeakNormalizedOperand {
  fn new(neg_log: &Array1<f64>) -> Self {
    let peak = min_or(neg_log, f64::INFINITY);
    if !peak.is_finite() {
      return Self {
        peak,
        relative_neg_log: Array1::from_elem(neg_log.len(), f64::INFINITY),
        plain: Array1::zeros(neg_log.len()),
      };
    }
    let relative_neg_log = neg_log.mapv(|nl| nl - peak);
    let plain = relative_neg_log.mapv(|r| (-r).exp());
    Self {
      peak,
      relative_neg_log,
      plain,
    }
  }

  fn tilted(&self, theta: f64, dx: f64) -> Array1<f64> {
    Zip::indexed(&self.relative_neg_log).map_collect(|i, &r| r - theta * (i.to_f64_exact() * dx))
  }
}

struct EndSlopes {
  left: f64,
  right: f64,
}

impl EndSlopes {
  fn of(neg_log: &Array1<f64>, margin: usize, dx: f64) -> Self {
    let n = neg_log.len();
    let width = margin.to_f64_exact() * dx;
    Self {
      left: (neg_log[margin] - neg_log[0]) / width,
      right: (neg_log[n - 1] - neg_log[n - 1 - margin]) / width,
    }
  }
}

struct TiltedConvolution {
  conv: Array1<f64>,
  floor: f64,
  offset: f64,
  theta: f64,
  dx: f64,
}

impl TiltedConvolution {
  fn new(operands: &ConvOperands, dx: f64, theta: f64) -> Result<Self, Report> {
    let za = operands.a.tilted(theta, dx);
    let zb = operands.b.tilted(theta, dx);
    let min_a = min_or(&za, f64::INFINITY);
    let min_b = min_or(&zb, f64::INFINITY);
    let pa = za.mapv(|z| (min_a - z).exp());
    let pb = zb.mapv(|z| (min_b - z).exp());
    let conv = convolve_fft(dx, &pa, &pb)?;
    let floor = max_or(&conv, 0.0) * CONV_TRUST_FRACTION;
    Ok(Self {
      conv,
      floor,
      offset: min_a + min_b,
      theta,
      dx,
    })
  }

  fn neg_log_at(&self, i: usize) -> Option<f64> {
    let c = self.conv[i];
    (c > self.floor).then(|| -c.ln() + self.offset + self.theta * (i.to_f64_exact() * self.dx))
  }
}

struct ResolvedNegLog {
  start: usize,
  neg_log: Array1<f64>,
}

struct ReconstructedConv {
  start_offset: usize,
  neg_log: Array1<f64>,
}
