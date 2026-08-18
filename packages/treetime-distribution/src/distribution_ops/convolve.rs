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

/// Fraction of the plain-space convolution peak below which FFT output is roundoff, not signal.
/// Values under this floor are discarded and the tail is reconstructed by fit instead. Matches v0
/// (`node_interpolator.py`: `fft_res > fft_res.max() * 1e-13`).
const CONV_TRUST_FRACTION: f64 = 1e-13;

/// Number of outermost trusted points used for the secant tail slope (v0 caps the margin at 3).
const CONV_TAIL_MARGIN: usize = 3;

pub fn distribution_convolution<Y: SupportsConvolution>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
    (Distribution::Formula(_), Distribution::Empty) | (Distribution::Empty, Distribution::Formula(_)) => {
      make_error!("Cannot convolve Formula with Empty: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Point(_)) | (Distribution::Point(_), Distribution::Formula(_)) => {
      make_error!("Cannot convolve Formula with Point: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Range(_)) | (Distribution::Range(_), Distribution::Formula(_)) => {
      make_error!("Cannot convolve Formula with Range: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Function(_)) | (Distribution::Function(_), Distribution::Formula(_)) => {
      make_error!("Cannot convolve Formula with Function: operation not implemented")
    },
    (Distribution::Formula(_), Distribution::Formula(_)) => {
      make_error!("Cannot convolve Formula with Formula: operation not implemented")
    },
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      guarded_empty_result("convolution", distribution_hard_domain(a), distribution_hard_domain(b))
    },
    (Distribution::Point(a), Distribution::Point(b)) => {
      Ok(convolution_point_point::<Y>(a, b)) //
    },
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      Ok(convolution_point_range::<Y>(a, b)) //
    },
    (Distribution::Range(a), Distribution::Range(b)) => {
      convolution_range_range::<Y>(a, b) //
    },
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      Ok(Distribution::Function(convolution_point_function::<Y>(a, b)?)) //
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      convolution_range_function::<Y>(a, b) //
    },
    (Distribution::Function(a), Distribution::Function(b)) => {
      convolution_function_function::<Y>(a, b) //
    },
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
  // Zero probability at the shoulders is the multiplicative identity, not a literal 0.0: it is
  // `0.0` under Plain but `+inf` under NegLog.
  let zero = Y::from_plain(0.0);

  if ulps_eq!(&peak_start, &peak_end, max_ulps = 10) {
    // Triangle case: equal-width ranges produce triangular shape
    let x = array![start, peak_start, end];
    let y = array![zero, peak_amplitude, zero];
    Distribution::function(x, y)
  } else {
    // Trapezoid case
    let x = array![start, peak_start, peak_end, end];
    let y = array![zero, peak_amplitude, peak_amplitude, zero];
    if has_uniform_spacing(&x) {
      // Trapezoid with uniform spacing
      Distribution::function(x, y)
    } else {
      // Trapezoid with non-uniform spacing: resample trapezoid to uniform grid
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

fn convolution_point_function<Y: SupportsConvolution>(
  p: &DistributionPoint<f64, Y>,
  f: &DistributionFunction<f64, Y>,
) -> Result<DistributionFunction<f64, Y>, Report> {
  let x_min = f.x_min() + p.t();
  let dx = f.dx();
  let y = f.y().mapv(|y| Y::multiply(y, p.amplitude()));
  DistributionFunction::from_start_dx_values(x_min, dx, y)
}

fn convolution_range_function<Y: SupportsConvolution>(
  r: &DistributionRange<f64, Y>,
  f: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  // DistributionFunction is guaranteed to be uniform
  let dx = f.dx();

  // split in a convolution with
  // - a point distribution (taking care of the shift + amplitude)
  // - an interval centered on zero and of a fixed width (taking care of the smoothing)
  let shift = f64::midpoint(r.start(), r.end());
  let amplitude = r.amplitude();
  let half_width = (r.end() - r.start()) / 2.0;

  let point_distr = DistributionPoint::new(shift, amplitude);
  let shifted_function = convolution_point_function::<Y>(&point_distr, f)?;

  // The box smoothing is an integral, so it runs in plain space. Peak-normalize to avoid NegLog
  // underflow, then restore the peak offset when converting the result back to the storage policy.
  let (plain, peak) = to_peak_normalized_plain::<Y>(shifted_function.y());
  let t = shifted_function.t();
  let mut smoothed = Array1::zeros(plain.len());
  // TODO: optimize by using cumulative sums
  for (i, &ti) in t.iter().enumerate() {
    let mask = t.mapv(|x| if (x - ti).abs() <= half_width { 1.0 } else { 0.0 });
    smoothed[i] = (&plain * &mask).sum() * dx;
  }
  let y = smoothed.mapv(|v| Y::from_neg_log(if v > 0.0 { -v.ln() + peak } else { f64::INFINITY }));

  DistributionFunction::from_start_dx_values(shifted_function.x_min(), shifted_function.dx(), y)
    .map(Distribution::Function)
}

/// Hard domain of a convolution operand for the empty-result guard: `None` when the operand carries
/// no mass (a legitimately empty convolution), otherwise the whole real line.
///
/// Convolution is a Minkowski sum, not an intersection: it translates an operand's mass across the
/// other operand's entire support, so two mass-bearing operands always jointly produce mass and are
/// never disjoint. Modelling a mass-bearing operand as unbounded on both sides (infinite bounds)
/// makes [`guarded_empty_result`] permit an empty result only when an operand itself is empty, and
/// flag any other empty convolution as the numerical collapse it is. The infinite bounds already
/// carry the unboundedness, so the per-side behavior is immaterial and left at the default.
fn conv_operand_domain(has_mass: bool) -> Option<HardDomain> {
  has_mass.then_some((
    (f64::NEG_INFINITY, f64::INFINITY),
    (BoundaryBehavior::Error, BoundaryBehavior::Error),
  ))
}

// Numerical routine: `a`/`b` operands, `t` grid, `n` trusted-point count follow the
// convolution/interpolation conventions used throughout this module.
#[allow(clippy::many_single_char_names)]
fn convolution_function_function<Y: SupportsConvolution>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  if a.is_empty() || b.is_empty() {
    return guarded_empty_result(
      "convolution",
      conv_operand_domain(!a.is_empty()),
      conv_operand_domain(!b.is_empty()),
    );
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
    return guarded_empty_result(
      "convolution",
      conv_operand_domain(!a.is_empty()),
      conv_operand_domain(!b.is_empty()),
    );
  }

  if a.len() == 1 && b.len() == 1 {
    return Ok(Distribution::point(
      a.x_min() + b.x_min(),
      Y::multiply(a.y()[0], b.y()[0]),
    ));
  }

  // A convolution is an integral, valid only in plain probability space. Convert each operand to
  // peak-normalized plain values (largest = 1) so NegLog's dynamic range cannot underflow; the peak
  // offsets are tracked in negative-log units and restored afterward.
  let (pa, peak_a) = to_peak_normalized_plain::<Y>(a.y());
  let (pb, peak_b) = to_peak_normalized_plain::<Y>(b.y());
  if !(peak_a.is_finite() && peak_b.is_finite()) {
    // A non-finite peak means that operand has no mass (all ordinates underflow to zero), so the
    // convolution is legitimately empty; the guard permits it via the massless operand.
    return guarded_empty_result(
      "convolution",
      conv_operand_domain(peak_a.is_finite()),
      conv_operand_domain(peak_b.is_finite()),
    );
  }

  let conv = convolve_fft(dx, &pa, &pb)?;

  // The FFT is trustworthy only in the bulk; below the roundoff floor the far tails are
  // reconstructed by log-linear extrapolation (v0 `NodeInterpolator.convolve_fft`).
  let Some(reconstructed) = reconstruct_neg_log_tails(&conv, dx)? else {
    // Both operands carry mass (finite peaks), so their convolution has mass: a collapse to empty
    // here is numerical, never structural. Both operands are present (unbounded for the sum), so the
    // guard reports the internal error rather than silently returning empty.
    return guarded_empty_result("convolution", conv_operand_domain(true), conv_operand_domain(true));
  };

  // The raw convolution starts at `a.x_min() + b.x_min()`; the reconstruction crops leading
  // untrusted samples, so the kept grid starts `start_offset` cells later.
  let x_min = a.x_min() + b.x_min() + (reconstructed.start_offset as f64) * dx;

  // Restore the peak offsets (negative-log units) and convert to the storage policy. `from_neg_log`
  // keeps NegLog's dynamic range and collapses sub-underflow values to zero under Plain.
  let offset = peak_a + peak_b;
  let y = reconstructed.neg_log.mapv(|v| Y::from_neg_log(v + offset));

  // A convolution whose trusted region collapses to a single cell is a point mass.
  if y.len() == 1 {
    return Ok(Distribution::point(x_min, y[0]));
  }

  let conv_distr = DistributionFunction::<f64, Y>::from_start_dx_values(x_min, dx, y)?;
  coarsen_convolution(conv_distr, dx_a.max(dx_b))
}

/// Coarsen the fine-grid convolution result back to the coarser operand spacing.
///
/// The FFT runs on the finest operand spacing; coarsening the result to the coarser operand spacing
/// keeps the point count near operand resolution. A trusted bulk narrower than half a coarse cell
/// would resample to a single point, which `Grid::from_range_dx` (and thus `resample_dx`) rejects.
/// In that degenerate case keep the fine-grid result unchanged: it is already a valid `>= 2`-point
/// Function, and returning it preserves the full shape instead of failing. This matches v0, whose
/// `convolve_fft` never coarsens and always returns the fine interpolation object.
pub(super) fn coarsen_convolution<Y: SupportsConvolution>(
  conv_distr: DistributionFunction<f64, Y>,
  coarse_dx: f64,
) -> Result<Distribution<Y>, Report> {
  // Point count of a spacing-`coarse_dx` grid over the result range, matching `Grid::from_range_dx`
  // (`round(range / dx) + 1`). The caller guarantees `y.len() >= 2`, so the range is positive.
  let range = conv_distr.x_max() - conv_distr.x_min();
  let coarse_points = (range / coarse_dx).round() as usize + 1;
  if coarse_points < 2 {
    return Ok(Distribution::Function(conv_distr));
  }
  Ok(Distribution::Function(conv_distr.resample_dx(coarse_dx)?))
}

/// Convert stored ordinates to peak-normalized plain probability for the convolution integral.
///
/// Returns the plain array (largest value 1.0, so no underflow) and the peak in negative-log units,
/// used to restore the offset after the FFT. Under a negative-log policy this is a
/// shift-and-exponentiate; under plain storage it reduces to dividing by the maximum.
fn to_peak_normalized_plain<Y: SupportsConvolution>(y: &Array1<f64>) -> (Array1<f64>, f64) {
  let neg_log = y.mapv(|v| Y::to_neg_log(v));
  let peak = min_or(&neg_log, f64::INFINITY);
  if !peak.is_finite() {
    return (Array1::zeros(y.len()), peak);
  }
  let plain = neg_log.mapv(|nl| (peak - nl).exp());
  (plain, peak)
}

/// Trusted-bulk reconstruction of a plain-space convolution in peak-relative negative-log space.
///
/// `neg_log` holds the reconstructed ordinates starting at grid index `start_offset` of the raw
/// convolution: the trusted bulk (`-ln conv` where `conv` sits above the roundoff floor) plus the
/// log-linearly extrapolated decaying tails. Non-decaying tail regions are cropped rather than
/// stored, so no `+inf` ordinate is ever produced (which would zero a whole interpolation cell).
struct ReconstructedConv {
  /// Index of the first kept sample within the raw convolution grid.
  start_offset: usize,
  /// Reconstructed negative-log ordinates over the kept range.
  neg_log: Array1<f64>,
}

/// Reconstruct the convolution result in peak-relative negative-log space, rebuilding the tails that
/// fall below the FFT roundoff floor by log-linear extrapolation.
///
/// Mirrors v0's `NodeInterpolator.convolve_fft`: crop to the trusted region (`conv > peak · 1e-13`),
/// then extend each side by the two-point secant slope over the outermost [`CONV_TAIL_MARGIN`]
/// trusted points, but only where the tail decays away from support (negative-log rising outward:
/// left slope negative toward smaller `t`, right slope positive toward larger `t`). A non-decaying
/// side is cropped at the bulk edge; the caller's boundary policy governs the domain beyond it.
/// Returns `None` when no trusted signal survives.
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

  // Slopes need only coordinate differences, so the relative grid origin is immaterial: adjacent
  // trusted points are `dx` apart and the secant spans `margin` cells.
  let left_slope = (trusted_y[margin] - trusted_y[0]) / (margin as f64 * dx);
  let right_slope = (trusted_y[n - 1] - trusted_y[n - 1 - margin]) / (margin as f64 * dx);

  // Extend a side only where the tail decays away from support and untrusted samples remain to
  // cover; otherwise crop at the bulk edge.
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
