use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::distribution::y_axis_policy::{Plain, SupportsConvolution};
use crate::make_error;
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, array};
use treetime_convolution::convolve;
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

#[cfg(test)]
mod tests {
  use super::distribution_convolution;
  use crate::distribution::distribution::DistributionPlain as Distribution;
  use ndarray::array;

  use pretty_assertions::assert_eq;

  #[test]
  fn test_convolution_empty() {
    let a: Distribution = Distribution::empty();
    let b: Distribution = Distribution::function(array![], array![]).unwrap();
    let actual: Distribution = distribution_convolution(&a, &b).unwrap();
    let expected: Distribution = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_point() {
    let a: Distribution = Distribution::point(2.0, 3.0);
    let b: Distribution = Distribution::point(5.0, 4.0);
    let actual: Distribution = distribution_convolution(&a, &b).unwrap();
    let expected: Distribution = Distribution::point(7.0, 12.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_triangle() {
    // Triangle case: equal-width ranges always produce triangles with uniform spacing
    let a = Distribution::range((2.0, 4.0), 3.0);
    let b = Distribution::range((6.0, 8.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      // start = 2.0 + 6.0 = 8.0
      // end = 4.0 + 8.0 = 12.0
      // peak_start = max(2.0+6.0, 4.0+6.0) = 10.0
      // peak_end = min(4.0+8.0, 2.0+8.0) = 10.0
      // Peak amplitude = 3.0 * 2.0 = 6.0
      let x = array![8.0, 10.0, 12.0];
      let y = array![0.0, 6.0, 0.0];
      Distribution::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_trapezoid_non_uniform() {
    let a = Distribution::range((2.0, 4.0), 3.0);
    let b = Distribution::range((6.0, 9.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      // Trapezoid: start=8.0, peak_start=10.0, peak_end=11.0, end=13.0
      // Non-uniform spacing (1, 2, 1.5, 2), so resampled to uniform grid
      // With dx=1.0 (smallest spacing), the uniform grid is:
      let x = array![8.0, 9.0, 10.0, 11.0, 12.0, 13.0];
      let y = array![0.0, 3.0, 6.0, 6.0, 3.0, 0.0];
      Distribution::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_trapezoid_uniform() {
    let a = Distribution::range((0.0, 2.0), 1.0);
    let b = Distribution::range((3.0, 7.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      // Trapezoid: start=3.0, peak_start=5.0, peak_end=7.0, end=9.0
      // Uniform spacing (dx=2.0), so no resampling needed
      let x = array![3.0, 5.0, 7.0, 9.0];
      let y = array![0.0, 2.0, 2.0, 0.0];
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
    let b_y = array![1.0, 2.0];
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
    let parent_y = array![0.2, 0.6, 0.2];
    let parent_dist = Distribution::function(parent_x, parent_y).unwrap();

    let branch_x = array![1.0, 1.5, 2.0];
    let branch_y = array![0.3, 0.4, 0.3];
    let branch_dist = Distribution::function(branch_x, branch_y).unwrap();

    let actual = distribution_convolution(&parent_dist, &branch_dist).unwrap();

    let expected_x = array![2011.0, 2011.5, 2012.0, 2012.5, 2013.0];
    let expected_y = array![0.03, 0.13, 0.18, 0.13, 0.03];
    let expected = Distribution::function(expected_x, expected_y).unwrap();

    assert_eq!(expected, actual);
  }
}
