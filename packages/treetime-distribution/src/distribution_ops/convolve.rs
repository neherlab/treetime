use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::distribution_ops::multiply::{HardDomain, distribution_hard_domain, guarded_empty_result};
use crate::policy::SupportsConvolution;
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, array, s};
use treetime_grid::BoundaryBehavior;
use treetime_ops::convolve_fft;
use treetime_utils::array::ndarray::{has_uniform_spacing, max_or, min_or};
use treetime_utils::make_error;

const CONV_TRUST_FRACTION: f64 = 1e-13;

const CONV_TAIL_MARGIN: usize = 3;

pub(crate) fn distribution_convolution_fine<Y: SupportsConvolution>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Function(a), Distribution::Function(b)) => convolution_function_function_fine::<Y>(a, b),
    _ => distribution_convolution(a, b),
  }
}

pub(crate) fn distribution_convolution<Y: SupportsConvolution>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Formula(_), _) | (_, Distribution::Formula(_)) => {
      make_error!("Cannot convolve {a} with {b}: operation not implemented")
    },
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      let a_domain = distribution_hard_domain(a);
      let b_domain = distribution_hard_domain(b);
      guarded_empty_result("convolution", a_domain, b_domain)
    },
    (Distribution::Point(a), Distribution::Point(b)) => Ok(convolution_point_point::<Y>(a, b)),
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      Ok(convolution_point_range::<Y>(a, b))
    },
    (Distribution::Range(a), Distribution::Range(b)) => convolution_range_range::<Y>(a, b),
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      Ok(Distribution::Function(convolution_point_function::<Y>(a, b)?))
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      convolution_range_function::<Y>(a, b)
    },
    (Distribution::Function(a), Distribution::Function(b)) => convolution_function_function::<Y>(a, b),
  }
}

fn convolution_point_point<Y: SupportsConvolution>(
  a: &DistributionPoint<f64, Y>,
  b: &DistributionPoint<f64, Y>,
) -> Distribution<Y> {
  let x = a.t() + b.t();
  let y = Y::multiply(a.amplitude(), b.amplitude());
  Distribution::point(x, y)
}

fn convolution_range_range<Y: SupportsConvolution>(
  a: &DistributionRange<f64, Y>,
  b: &DistributionRange<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let start = a.start() + b.start();
  let end = a.end() + b.end();

  let peak_start = f64::max(a.start() + b.start(), a.end() + b.start());
  let peak_end = f64::min(a.end() + b.end(), a.start() + b.end());

  let peak_amplitude = Y::multiply(a.amplitude(), b.amplitude());
  let zero = Y::from_plain(0.0);

  if ulps_eq!(&peak_start, &peak_end, max_ulps = 10) {
    let x = array![start, peak_start, end];
    let y = array![zero, peak_amplitude, zero];
    Distribution::function(x, y)
  } else {
    let x = array![start, peak_start, peak_end, end];
    let y = array![zero, peak_amplitude, peak_amplitude, zero];
    if has_uniform_spacing(&x) {
      Distribution::function(x, y)
    } else {
      DistributionFunction::from_arrays_nonuniform(&x, &y).map(Distribution::Function)
    }
  }
}

fn convolution_point_range<Y: SupportsConvolution>(
  p: &DistributionPoint<f64, Y>,
  r: &DistributionRange<f64, Y>,
) -> Distribution<Y> {
  let begin = r.start() + p.t();
  let end = r.end() + p.t();
  let amplitude = Y::multiply(p.amplitude(), r.amplitude());
  Distribution::range((begin, end), amplitude)
}

fn convolution_range_function<Y: SupportsConvolution>(
  r: &DistributionRange<f64, Y>,
  f: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let dx = f.dx();

  let shift = f64::midpoint(r.start(), r.end());
  let amplitude = r.amplitude();
  let half_width = (r.end() - r.start()) / 2.0;

  let point_distr = DistributionPoint::new(shift, amplitude);
  let shifted_function = convolution_point_function::<Y>(&point_distr, f)?;

  let (plain, peak) = to_peak_normalized_plain::<Y>(shifted_function.y());
  let t = shifted_function.t();
  let mut smoothed = Array1::zeros(plain.len());
  for (i, &ti) in t.iter().enumerate() {
    let mask = t.mapv(|x| if (x - ti).abs() <= half_width { 1.0 } else { 0.0 });
    smoothed[i] = (&plain * &mask).sum() * dx;
  }
  let y = smoothed.mapv(|v| Y::from_neg_log(if v > 0.0 { -v.ln() + peak } else { f64::INFINITY }));

  DistributionFunction::from_start_dx_values(shifted_function.x_min(), shifted_function.dx(), y)
    .map(Distribution::Function)
}

fn convolution_point_function<Y: SupportsConvolution>(
  p: &DistributionPoint<f64, Y>,
  f: &DistributionFunction<f64, Y>,
) -> Result<DistributionFunction<f64, Y>, Report> {
  let x_min = f.x_min() + p.t();
  let dx = f.dx();
  let y = f.y().mapv(|y| Y::multiply(y, p.amplitude()));
  DistributionFunction::from_start_dx_values(x_min, dx, y)
}

fn convolution_function_function<Y: SupportsConvolution>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  let coarse_dx = a.dx().max(b.dx());
  match convolution_function_function_fine(a, b)? {
    Distribution::Function(conv_distr) => coarsen_convolution(conv_distr, coarse_dx),
    other => Ok(other),
  }
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
fn convolution_function_function_fine<Y: SupportsConvolution>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
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

  let a = a.resample_dx(dx)?;
  let b = b.resample_dx(dx)?;
  if a.is_empty() || b.is_empty() {
    let a_domain = conv_operand_domain(!a.is_empty());
    let b_domain = conv_operand_domain(!b.is_empty());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  if a.len() == 1 && b.len() == 1 {
    return Ok(Distribution::point(
      a.x_min() + b.x_min(),
      Y::multiply(a.y()[0], b.y()[0]),
    ));
  }

  let (pa, peak_a) = to_peak_normalized_plain::<Y>(a.y());
  let (pb, peak_b) = to_peak_normalized_plain::<Y>(b.y());
  if !(peak_a.is_finite() && peak_b.is_finite()) {
    let a_domain = conv_operand_domain(peak_a.is_finite());
    let b_domain = conv_operand_domain(peak_b.is_finite());
    return guarded_empty_result("convolution", a_domain, b_domain);
  }

  let conv = convolve_fft(dx, &pa, &pb)?;

  let Some(reconstructed) = reconstruct_neg_log_tails(&conv, dx)? else {
    let a_domain = conv_operand_domain(true);
    let b_domain = conv_operand_domain(true);
    return guarded_empty_result("convolution", a_domain, b_domain);
  };

  let x_min = a.x_min() + b.x_min() + (reconstructed.start_offset as f64) * dx;

  let offset = peak_a + peak_b;
  let y = reconstructed.neg_log.mapv(|v| Y::from_neg_log(v + offset));

  if y.len() == 1 {
    return Ok(Distribution::point(x_min, y[0]));
  }

  let conv_distr = DistributionFunction::<f64, Y>::from_start_dx_values(x_min, dx, y)?;
  Ok(Distribution::Function(conv_distr))
}

fn conv_operand_domain(has_mass: bool) -> Option<HardDomain> {
  has_mass.then_some((
    (f64::NEG_INFINITY, f64::INFINITY),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
  ))
}

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub(super) fn coarsen_convolution<Y: SupportsConvolution>(
  conv_distr: DistributionFunction<f64, Y>,
  coarse_dx: f64,
) -> Result<Distribution<Y>, Report> {
  let range = conv_distr.x_max() - conv_distr.x_min();
  let coarse_points = (range / coarse_dx).round() as usize + 1;
  if coarse_points < 2 {
    return Ok(Distribution::Function(conv_distr));
  }
  Ok(Distribution::Function(conv_distr.resample_dx(coarse_dx)?))
}

fn to_peak_normalized_plain<Y: SupportsConvolution>(y: &Array1<f64>) -> (Array1<f64>, f64) {
  let neg_log = y.mapv(|v| Y::to_neg_log(v));
  let peak = min_or(&neg_log, f64::INFINITY);
  if !peak.is_finite() {
    return (Array1::zeros(y.len()), peak);
  }
  let plain = neg_log.mapv(|nl| (peak - nl).exp());
  (plain, peak)
}

#[allow(
  clippy::as_conversions,
  clippy::integer_division,
  reason = "count/index numeric cast is exact for the domain range; integer division is the intended floor division"
)]
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

  let left_slope = (trusted_y[margin] - trusted_y[0]) / (margin as f64 * dx);
  let right_slope = (trusted_y[n - 1] - trusted_y[n - 1 - margin]) / (margin as f64 * dx);

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
    neg_log[i - left_start] = trusted_y[0] + left_slope * ((i as f64 - first as f64) * dx);
  }
  for i in (last + 1)..=right_end {
    neg_log[i - left_start] = trusted_y[n - 1] + right_slope * ((i as f64 - last as f64) * dx);
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
