use crate::Distribution;
use crate::distribution_core::function::DistributionFunction;
use crate::distribution_core::point::DistributionPoint;
use crate::distribution_core::range::DistributionRange;
use crate::policy::{NegLog, Plain, SupportsConvolution};
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, array, s};
use treetime_ops::convolve;
use treetime_utils::array::ndarray::{has_uniform_spacing, max_or, min_or};
use treetime_utils::make_error;

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
      Ok(Distribution::Empty) //
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

  if ulps_eq!(&peak_start, &peak_end, max_ulps = 10) {
    // Triangle case: equal-width ranges produce triangular shape
    let x = array![start, peak_start, end];
    let y = array![0.0, peak_amplitude, 0.0];
    Distribution::function(x, y)
  } else {
    // Trapezoid case
    let x = array![start, peak_start, peak_end, end];
    let y = array![0.0, peak_amplitude, peak_amplitude, 0.0];
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

  // Convolution with a range centered on zero and of given width
  let mut y_out = Array1::zeros(shifted_function.y().len());

  // TODO: optimize by using cumulative sums
  for (i, &ti) in shifted_function.t().iter().enumerate() {
    let mask = shifted_function.t().mapv(|x| (x - ti).abs() <= half_width);
    let filtered_y = shifted_function.y() * &mask.mapv(|x| if x { 1.0 } else { 0.0 });
    y_out[i] = filtered_y.sum() * dx;
  }

  DistributionFunction::from_start_dx_values(shifted_function.x_min(), shifted_function.dx(), y_out)
    .map(Distribution::Function)
}

fn convolution_function_function<Y: SupportsConvolution>(
  a: &DistributionFunction<f64, Y>,
  b: &DistributionFunction<f64, Y>,
) -> Result<Distribution<Y>, Report> {
  // Check for degenerate cases
  if a.is_empty() || b.is_empty() {
    return Ok(Distribution::empty());
  }

  let dx_a = a.dx();
  let dx_b = b.dx();
  let dx = dx_a.min(dx_b);

  if !(dx.is_finite() && dx > 0.0) {
    return make_error!("Invalid grid spacing detected during convolution: {dx}");
  }

  let a_resampled = a.resample_dx(dx)?;
  let b_resampled = b.resample_dx(dx)?;

  if a_resampled.is_empty() || b_resampled.is_empty() {
    return Ok(Distribution::empty());
  }

  if a_resampled.len() == 1 && b_resampled.len() == 1 {
    return Ok(Distribution::point(
      a_resampled.x_min() + b_resampled.x_min(),
      Y::multiply(a_resampled.y()[0], b_resampled.y()[0]),
    ));
  }

  let output_len = a_resampled.len() + b_resampled.len() - 1;
  let output_grid_start = a_resampled.x_min() + b_resampled.x_min();

  let conv_result = convolve(dx, a_resampled.y(), b_resampled.y())?;
  if conv_result.len() != output_len {
    return make_error!(
      "Convolution result length mismatch: expected {}, got {}",
      output_len,
      conv_result.len()
    );
  }
  let conv_distr = DistributionFunction::<f64, Plain>::from_start_dx_values(output_grid_start, dx, conv_result)?;

  let coarse_dx = f64::max(dx_a, dx_b);
  let final_distr = conv_distr.resample_dx(coarse_dx)?;

  if final_distr.len() < 2 {
    return make_error!("Final distribution after convolution has less than two points");
  }

  let final_distr = DistributionFunction::<f64, Y>::from_start_dx_values(
    final_distr.x_min(),
    final_distr.dx(),
    final_distr.y().clone(),
  )?;
  Ok(Distribution::Function(final_distr))
}

// ---------------------------------------------------------------------------------------------
// DRAFT: log-space (NegLog) convolution round-trip. Not yet wired into the pipeline. Scope is the
// Function * Function case, which carries the FFT and the tail-reconstruction science. Point/range
// convolutions need no FFT and are deferred.
// ---------------------------------------------------------------------------------------------

/// Fraction of the plain-space convolution peak below which FFT output is roundoff, not signal.
/// Values under this floor are discarded and the tail is reconstructed by fit instead. Matches v0
/// (`node_interpolator.py`: `fft_res > fft_res.max() * 1e-13`).
const CONV_TRUST_FRACTION: f64 = 1e-13;

/// Number of outermost trusted points used for the secant tail slope (v0 caps the margin at 3).
const CONV_TAIL_MARGIN: usize = 3;

/// Convolve two negative-log distributions.
///
/// A convolution is an integral, not a pointwise operation, so it cannot be carried out in
/// negative-log space. Each operand is moved to plain probability (peak-subtracted so the largest
/// value is 1.0 and nothing underflows), convolved by FFT, then returned to negative-log space with
/// the two peak offsets restored. The FFT is trustworthy only in the bulk: values below
/// [`CONV_TRUST_FRACTION`] of the peak are roundoff. The decisive far-tail values (the node-603
/// case is ~1e-12 of peak, well under the roundoff floor) are therefore **reconstructed by
/// log-linear extrapolation** from the outermost trusted points rather than read out of the FFT,
/// following v0's `NodeInterpolator.convolve_fft`.
///
/// DRAFT: implements Function * Function only.
pub fn distribution_convolution_neglog(
  a: &Distribution<NegLog>,
  b: &Distribution<NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => Ok(Distribution::Empty),
    (Distribution::Function(a), Distribution::Function(b)) => convolution_function_function_neglog(a, b),
    _ => make_error!("Log-space convolution is drafted for Function * Function only"),
  }
}

// Numerical routine: `a`/`b` operands, `t` grid, `y` ordinates, `n` trusted-point count follow the
// convolution/interpolation conventions used throughout this module.
#[allow(clippy::many_single_char_names)]
fn convolution_function_function_neglog(
  a: &DistributionFunction<f64, NegLog>,
  b: &DistributionFunction<f64, NegLog>,
) -> Result<Distribution<NegLog>, Report> {
  if a.is_empty() || b.is_empty() {
    return Ok(Distribution::empty());
  }

  let dx = a.dx().min(b.dx());
  if !(dx.is_finite() && dx > 0.0) {
    return make_error!("Invalid grid spacing detected during log-space convolution: {dx}");
  }

  // Resample onto the finer common spacing. Still negative-log, so this is linear interpolation of
  // -ln(p): the intended piecewise-log-linear representation.
  let a = a.resample_dx(dx)?;
  let b = b.resample_dx(dx)?;
  if a.is_empty() || b.is_empty() {
    return Ok(Distribution::empty());
  }

  // Peak-subtract each operand, then exponentiate to plain probability with peak 1.0. The shifts
  // are restored afterward as an additive constant in negative-log space.
  let shift_a = min_or(a.y(), f64::INFINITY);
  let shift_b = min_or(b.y(), f64::INFINITY);
  if !(shift_a.is_finite() && shift_b.is_finite()) {
    return Ok(Distribution::empty());
  }
  let pa = a.y().mapv(|y| (shift_a - y).exp());
  let pb = b.y().mapv(|y| (shift_b - y).exp());

  // Integral in plain space. `convolve` already applies dx (Riemann normalization).
  let conv = convolve(dx, &pa, &pb)?;
  let x_min = a.x_min() + b.x_min();
  let t = Array1::from_shape_fn(conv.len(), |i| x_min + (i as f64) * dx);
  let offset = shift_a + shift_b;

  // Trusted bulk: the contiguous region where the plain convolution exceeds the roundoff floor.
  let peak = max_or(&conv, 0.0);
  if peak <= 0.0 {
    return Ok(Distribution::empty());
  }
  let floor = peak * CONV_TRUST_FRACTION;
  let (Some(first), Some(last)) = (
    conv.iter().position(|&c| c > floor),
    conv.iter().rposition(|&c| c > floor),
  ) else {
    return Ok(Distribution::empty());
  };

  let trusted_t = t.slice(s![first..=last]).to_owned();
  let trusted_y = conv.slice(s![first..=last]).mapv(|c| -c.ln() + offset);
  let n = trusted_y.len();
  let margin = CONV_TAIL_MARGIN.min(n / 3);
  if margin < 1 {
    return make_error!("Log-space convolution left too few trusted points to reconstruct tails");
  }

  // v0-style secant slopes over the outermost `margin` trusted points. A tail is extended only
  // where it decays away from support: the neg-log ordinate rises outward (left slope < 0 toward
  // smaller t, right slope > 0 toward larger t). Where it does not decay, the region carries no
  // mass (+inf neg-log = zero probability).
  let left_slope = (trusted_y[margin] - trusted_y[0]) / (trusted_t[margin] - trusted_t[0]);
  let right_slope = (trusted_y[n - 1] - trusted_y[n - 1 - margin]) / (trusted_t[n - 1] - trusted_t[n - 1 - margin]);

  let mut y = Array1::from_elem(conv.len(), f64::INFINITY);
  for i in first..=last {
    y[i] = -conv[i].ln() + offset;
  }
  if left_slope < 0.0 {
    for i in 0..first {
      y[i] = trusted_y[0] + left_slope * (t[i] - trusted_t[0]);
    }
  }
  if right_slope > 0.0 {
    for i in (last + 1)..conv.len() {
      y[i] = trusted_y[n - 1] + right_slope * (t[i] - trusted_t[n - 1]);
    }
  }

  // TODO(Part B): carry the fitted slopes as soft `BoundaryBehavior::Linear` tails so the law
  // continues beyond the grid. Until then the reconstructed tail lives on the grid itself.
  let result = DistributionFunction::<f64, NegLog>::from_start_dx_values(x_min, dx, y)?;
  Ok(Distribution::Function(result).normalize())
}
