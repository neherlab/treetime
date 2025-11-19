use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::make_error;
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, array};
use treetime_convolution::convolve;
use treetime_utils::ndarray::ndarray_uniform_grid;

pub fn distribution_convolution(a: &Distribution, b: &Distribution) -> Result<Distribution, Report> {
  match (a, b) {
    (Distribution::Empty, _) | (_, Distribution::Empty) => {
      Ok(Distribution::Empty) //
    },
    (Distribution::Point(a), Distribution::Point(b)) => {
      Ok(convolution_point_point(a, b)) //
    },
    (Distribution::Point(a), Distribution::Range(b)) | (Distribution::Range(b), Distribution::Point(a)) => {
      Ok(convolution_point_range(a, b)) //
    },
    (Distribution::Range(a), Distribution::Range(b)) => {
      convolution_range_range(a, b) //
    },
    (Distribution::Point(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Point(a)) => {
      Ok(Distribution::Function(convolution_point_function(a, b)?))
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      convolution_range_function(a, b) //
    },
    (Distribution::Function(a), Distribution::Function(b)) => {
      convolution_function_function(a, b) //
    },
  }
}

fn convolution_point_point(a: &DistributionPoint<f64>, b: &DistributionPoint<f64>) -> Distribution {
  let x = a.t() + b.t();
  let y = a.amplitude() * b.amplitude();
  Distribution::point(x, y)
}

fn convolution_range_range(a: &DistributionRange<f64>, b: &DistributionRange<f64>) -> Result<Distribution, Report> {
  let start = a.start() + b.start();
  let end = a.end() + b.end();

  let peak_start = f64::max(a.start() + b.start(), a.end() + b.start());
  let peak_end = f64::min(a.end() + b.end(), a.start() + b.end());

  let peak_amplitude = a.amplitude() * b.amplitude();

  if ulps_eq!(&peak_start, &peak_end, max_ulps = 10) {
    let x = array![start, peak_start, end];
    let y = array![0.0, peak_amplitude, 0.0];
    Distribution::function(x, y)
  } else {
    let x = array![start, peak_start, peak_end, end];
    let y = array![0.0, peak_amplitude, peak_amplitude, 0.0];
    Distribution::function(x, y)
  }
}

fn convolution_point_range(p: &DistributionPoint<f64>, r: &DistributionRange<f64>) -> Distribution {
  let begin = r.start() + p.t();
  let end = r.end() + p.t();
  let amplitude = p.amplitude() * r.amplitude();
  Distribution::range((begin, end), amplitude)
}

fn convolution_point_function(
  p: &DistributionPoint<f64>,
  f: &DistributionFunction<f64>,
) -> Result<DistributionFunction<f64>, Report> {
  let x_min = f.x_min() + p.t();
  let dx = f.dx();
  let y = f.y().map(|y| y * p.amplitude());
  DistributionFunction::from_start_dx_values(x_min, dx, y)
}

fn convolution_range_function(
  r: &DistributionRange<f64>,
  f: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  // DistributionFunction is guaranteed to be uniform
  let dx = f.dx();

  // split in a convolution with
  // - a point distribution (taking care of the shift + amplitude)
  // - an interval centered on zero and of a fixed width (taking care of the smoothing)
  let shift = f64::midpoint(r.start(), r.end());
  let amplitude = r.amplitude();
  let half_width = (r.end() - r.start()) / 2.0;

  let point_distr = DistributionPoint::new(shift, amplitude);
  let shifted_function = convolution_point_function(&point_distr, f).unwrap();

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

fn convolution_function_function(
  a: &DistributionFunction<f64>,
  b: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  // Check for degenerate cases
  if a.is_empty() || b.is_empty() {
    return Ok(Distribution::empty());
  }

  let dx_a = a.dx();
  let dx_b = b.dx();
  // Use max instead of min to prevent exponentially growing grid sizes
  let dx = dx_a.min(dx_b);

  if !(dx.is_finite() && dx > 0.0) {
    return make_error!("Invalid grid spacing detected during convolution: {dx}");
  }

  let (a_values, a_min) = resample_distribution(a, dx)?;
  let (b_values, b_min) = resample_distribution(b, dx)?;

  // Handle empty distributions after resampling
  if a_values.len() == 0 || b_values.len() == 0 {
    return Ok(Distribution::empty());
  }

  // For single-point distributions, treat as delta functions
  if a_values.len() == 1 && b_values.len() == 1 {
    return Ok(Distribution::point(a_min + b_min, a_values[0] * b_values[0]));
  }

  let output_len = a_values.len() + b_values.len() - 1;
  let output_grid_start = a_min + b_min;

  // Perform convolution with unified grids
  let conv_result = convolve(dx, &a_values, &b_values)?;
  // check that the convolution result has the expected length
  if conv_result.len() != output_len {
    return make_error!(
      "Convolution result length mismatch: expected {}, got {}",
      output_len,
      conv_result.len()
    );
  }
  let conv_distr = DistributionFunction::from_start_dx_values(output_grid_start, dx, conv_result)?;

  // resample distribution to coarser grid
  let coarse_dx = dx_a.max(dx_b);
  let (final_values, coarse_grid_min) = resample_distribution(&conv_distr, coarse_dx)?;
  let coarse_grid = ndarray_uniform_grid(coarse_grid_min, coarse_dx, final_values.len());

  // sanity check: more than one point
  if coarse_grid.len() < 2 {
    return make_error!("Final distribution after convolution has less than two points");
  }

  // return final distribution
  Distribution::function(coarse_grid, final_values)
}

fn resample_distribution(func: &DistributionFunction<f64>, dx: f64) -> Result<(Array1<f64>, f64), Report> {
  if func.is_empty() {
    return make_error!("DistributionFunction is empty");
  }
  let min = func.x_min();
  let max = func.x_max();

  if max < min {
    return make_error!("DistributionFunction grid must be sorted ascending");
  }

  if (max - min).abs() < f64::EPSILON {
    // Single point case
    let value = func.interp(min)?;
    return Ok((array![value], min));
  }

  let original_span = max - min;

  // Calculate number of steps needed to cover the span with the given dx
  let steps = (original_span / dx).ceil() as usize;
  let len = steps + 1;

  if len < 2 {
    return make_error!("Unable to resample distribution: computed length {len} is too small");
  }

  // Prevent excessive memory allocation
  if len > 1_000_000 {
    return make_error!("Resampling would require {len} points, which exceeds safety limit");
  }

  let mut values = Array1::zeros(len);
  for i in 0..len {
    let position = min + dx * (i as f64);
    // Clamp position to valid range to avoid extrapolation errors
    let clamped_position = position.min(max).max(min);
    values[i] = func.interp(clamped_position)?;
  }

  Ok((values, min))
}

#[cfg(test)]
mod tests {
  use super::*;

  use pretty_assertions::assert_eq;

  #[test]
  fn test_convolution_empty() {
    let a = Distribution::empty();
    let b = Distribution::function(array![], array![]).unwrap();
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_point() {
    let a = Distribution::point(2.0, 3.0);
    let b = Distribution::point(5.0, 4.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::point(7.0, 12.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_same_length() {
    let a = Distribution::range((2.0, 4.0), 3.0);
    let b = Distribution::range((6.0, 8.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![8.0, 10.0, 12.0];
      let y = array![0.0, 6.0, 0.0];
      Distribution::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_different_length() {
    let a = Distribution::range((2.0, 4.0), 3.0);
    let b = Distribution::range((6.0, 9.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![8.0, 10.0, 11.0, 13.0];
      let y = array![0.0, 6.0, 6.0, 0.0];
      Distribution::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_function() {
    let p = Distribution::point(3.0, 2.0);

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let f = Distribution::function(x, y).unwrap();
    let actual = distribution_convolution(&p, &f).unwrap();

    let x = array![3.0, 4.0, 5.0, 6.0, 7.0];
    let y = array![2.0, 4.0, 6.0, 8.0, 10.0];
    let expected = Distribution::function(x, y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_function() {
    let r = Distribution::range((2.0, 6.0), 0.5);

    let x = array![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
    let y = array![0.0, 1.0, 0.0, 2.0, 1.0, 0.0];
    let f = Distribution::function(x, y).unwrap();
    let actual = distribution_convolution(&r, &f).unwrap();

    let x = array![4.0, 6.0, 8.0, 10.0, 12.0, 14.0];
    let y = array![1.0, 1.0, 3.0, 3.0, 3.0, 1.0];
    let expected = Distribution::function(x, y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_range() {
    let p = Distribution::point(3.0, 2.0);
    let r = Distribution::range((1.0, 4.0), 1.5);
    let actual = distribution_convolution(&p, &r).unwrap();
    let expected = Distribution::range((4.0, 7.0), 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_point() {
    let r = Distribution::range((1.0, 4.0), 1.5);
    let p = Distribution::point(3.0, 2.0);
    let actual = distribution_convolution(&r, &p).unwrap();
    let expected = Distribution::range((4.0, 7.0), 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_basic() {
    // Simple test with aligned grids
    let a_x = array![0.0, 1.0, 2.0];
    let a_y = array![1.0, 2.0, 1.0];
    let a = Distribution::function(a_x, a_y).unwrap();

    let b_x = array![0.0, 1.0];
    let b_y = array![1.0, 2.0]; // Make values non-uniform to force Function type
    let b = Distribution::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();

    let expected_x = array![0.0, 1.0, 2.0, 3.0];
    let expected_y = array![1.0, 4.0, 5.0, 2.0];
    let expected = Distribution::function(expected_x, expected_y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_single_points() {
    let a_x = array![2.0];
    let a_y = array![3.0];
    let a = Distribution::function(a_x, a_y).unwrap();

    let b_x = array![5.0];
    let b_y = array![4.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::point(7.0, 12.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_empty() {
    let a = Distribution::function(array![], array![]).unwrap();
    let b_x = array![1.0, 2.0];
    let b_y = array![1.0, 1.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_different_spacing() {
    // Test with different grid spacings to ensure proper resampling
    let a_x = array![0.0, 0.5, 1.0];
    let a_y = array![1.0, 2.0, 1.0];
    let a = Distribution::function(a_x, a_y).unwrap();

    let b_x = array![0.0, 1.0, 2.0];
    let b_y = array![1.0, 1.0, 1.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();

    let expected_x = array![0.0, 1.0, 2.0, 3.0];
    let expected_y = array![0.5, 2.0, 2.0, 0.5];
    let expected = Distribution::function(expected_x, expected_y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_zero_width() {
    // Test point distribution convolved with function
    let a = Distribution::point(5.0, 1.0);

    let b_x = array![1.0, 2.0];
    let b_y = array![2.0, 3.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();

    let expected_x = array![6.0, 7.0];
    let expected_y = array![2.0, 3.0];
    let expected = Distribution::function(expected_x, expected_y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_backward_pass_temporal_direction() {
    // Test the specific use case in backward pass: parent_time = child_time + (-branch_length)
    let child_time_dist = Distribution::point(2013.0, 1.0);
    let branch_length_dist = Distribution::point(2.5, 1.0);

    // In backward pass, we negate the branch distribution
    let negated_branch = branch_length_dist.negate();
    let actual = distribution_convolution(&child_time_dist, &negated_branch).unwrap();

    let expected = Distribution::point(2010.5, 1.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_forward_pass_temporal_direction() {
    // Test the forward pass: child_time = parent_time + branch_length
    let parent_time_dist = Distribution::point(2010.0, 1.0);
    let branch_length_dist = Distribution::point(1.5, 1.0);

    let actual = distribution_convolution(&parent_time_dist, &branch_length_dist).unwrap();

    let expected = Distribution::point(2011.5, 1.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_with_uncertainty() {
    // Test convolution with uncertainty distributions (functions)
    let parent_x = array![2010.0, 2010.5, 2011.0];
    let parent_y = array![0.2, 0.6, 0.2]; // Uncertainty around 2010.5
    let parent_dist = Distribution::function(parent_x, parent_y).unwrap();

    let branch_x = array![1.0, 1.5, 2.0];
    let branch_y = array![0.3, 0.4, 0.3]; // Uncertainty around 1.5
    let branch_dist = Distribution::function(branch_x, branch_y).unwrap();

    let actual = distribution_convolution(&parent_dist, &branch_dist).unwrap();

    let expected_x = array![2011.0, 2011.5, 2012.0, 2012.5, 2013.0];
    let expected_y = array![0.03, 0.13, 0.18, 0.13, 0.03];
    let expected = Distribution::function(expected_x, expected_y).unwrap();

    assert_eq!(expected, actual);
  }
}
