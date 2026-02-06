use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::distribution::y_axis_policy::{Plain, SupportsConvolution};
use crate::make_error;
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, array};
use treetime_ops::convolve;
use treetime_utils::ndarray::has_uniform_spacing;

pub fn distribution_convolution<Y: SupportsConvolution>(
  a: &Distribution<Y>,
  b: &Distribution<Y>,
) -> Result<Distribution<Y>, Report> {
  match (a, b) {
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
    (Distribution::Formula(_), _) | (_, Distribution::Formula(_)) => {
      panic!("Convolution not implemented for Formula distributions") //
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
  let t_out = shifted_function.t();
  let mut y_out = Array1::zeros(shifted_function.y().len());

  // TODO: optimize by using cumulative sums
  for (i, &ti) in shifted_function.t().iter().enumerate() {
    let mask = shifted_function.t().mapv(|x| (x - ti).abs() <= half_width);
    let filtered_y = shifted_function.y() * &mask.mapv(|x| if x { 1.0 } else { 0.0 });
    y_out[i] = filtered_y.sum() * dx;
  }

  Distribution::function(t_out, y_out)
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

