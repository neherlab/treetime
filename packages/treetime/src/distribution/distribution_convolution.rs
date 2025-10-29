use crate::distribution::distribution::Distribution;
use crate::distribution::distribution_function::DistributionFunction;
use crate::distribution::distribution_point::DistributionPoint;
use crate::distribution::distribution_range::DistributionRange;
use crate::make_error;
use crate::make_report;
use approx::ulps_eq;
use eyre::Report;
use ndarray::{Array1, array, s};
use std::cmp::max;
use treetime_convolution::convolve;
use treetime_utils::ndarray::{is_uniform_grid, ndarray_pad_zeros_right, ndarray_uniform_grid};

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
      Ok(convolution_point_function(a, b)?) //
    },
    (Distribution::Range(a), Distribution::Function(b)) | (Distribution::Function(b), Distribution::Range(a)) => {
      Ok(convolution_range_function(a, b)) //
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
) -> Result<Distribution, Report> {
  let t = f.t().map(|t| t + p.t());
  let y = f.y().map(|y| y * p.amplitude());
  Distribution::function(t, y)
}

fn convolution_range_function(r: &DistributionRange<f64>, f: &DistributionFunction<f64>) -> Distribution {
  let t_out = f.t().clone();
  let mut y_out = Array1::zeros(f.y().len());

  for (i, &ti) in f.t().iter().enumerate() {
    let mask = f.t().mapv(|x| (x >= ti - r.end()) && (x <= ti - r.start()));
    let filtered_y = f.y() * &mask.mapv(|x| if x { 1.0 } else { 0.0 });
    y_out[i] = r.amplitude() * filtered_y.sum();
  }

  Distribution::function(t_out, y_out).unwrap()
}

fn convolution_function_function(
  a: &DistributionFunction<f64>,
  b: &DistributionFunction<f64>,
) -> Result<Distribution, Report> {
  // Check for degenerate cases
  if a.t().is_empty() || b.t().is_empty() {
    return Ok(Distribution::empty());
  }

  let dx_a = compute_uniform_spacing(a.t())?;
  let dx_b = compute_uniform_spacing(b.t())?;
  // Use max instead of min to prevent exponentially growing grid sizes
  let dx = dx_a.max(dx_b);

  if !(dx.is_finite() && dx > 0.0) {
    return make_error!("Invalid grid spacing detected during convolution: {dx}");
  }

  let (a_values, a_min) = resample_distribution(a, dx)?;
  let (b_values, b_min) = resample_distribution(b, dx)?;

  if a_values.len() == 0 || b_values.len() == 0 {
    return Ok(Distribution::empty());
  }

  // For single-point distributions, treat as delta functions
  if a_values.len() == 1 && b_values.len() == 1 {
    return Ok(Distribution::point(a_min + b_min, a_values[0] * b_values[0]));
  }

  let output_len = a_values.len() + b_values.len() - 1;
  let output_grid_local = ndarray_uniform_grid(0.0, dx, output_len);

  let max_len = max(a_values.len(), b_values.len());
  let unified_grid = ndarray_uniform_grid(0.0, dx, max_len);

  // Resample both distributions onto the unified grid
  let a_unified = ndarray_pad_zeros_right(&a_values, max_len);
  let b_unified = ndarray_pad_zeros_right(&b_values, max_len);

  // Perform convolution with unified grids
  let unified_output_grid = ndarray_uniform_grid(0.0, dx, output_len);
  let unified_conv_result = convolve(&unified_grid, &a_unified, &b_unified, &unified_output_grid)?;

  // Extract only the meaningful portion of the result (first output_len elements)
  let conv_result = unified_conv_result.slice(s![..output_len]).to_owned();

  // Shift the output grid to the correct position
  let shift = a_min + b_min;
  let output_grid = output_grid_local.mapv(|x| x + shift);

  Distribution::function(output_grid, conv_result)
}

fn compute_uniform_spacing(grid: &Array1<f64>) -> Result<f64, Report> {
  if grid.len() < 2 {
    return make_error!("Function distributions require at least two grid points");
  }

  if !is_uniform_grid(grid) {
    return make_error!("Function distributions must use uniform grids for convolution");
  }

  let spacing = grid[1] - grid[0];
  if !(spacing.is_finite() && spacing > 0.0) {
    return make_error!("Encountered non-positive grid spacing: {spacing}");
  }

  Ok(spacing)
}

fn resample_distribution(func: &DistributionFunction<f64>, dx: f64) -> Result<(Array1<f64>, f64), Report> {
  let t = func.t();
  let min = *t.first().ok_or_else(|| make_report!("DistributionFunction is empty"))?;
  let max = *t.last().unwrap();

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
    let r = Distribution::range((1.0, 3.0), 2.0);

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![0.0, 1.0, 0.0, 2.0, 1.0, 0.0];
    let f = Distribution::function(x, y).unwrap();
    let actual = distribution_convolution(&r, &f).unwrap();

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![0.0, 0.0, 2.0, 2.0, 6.0, 6.0];
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

    let result = distribution_convolution(&a, &b).unwrap();

    match result {
      Distribution::Function(f) => {
        assert_eq!(f.t().len(), 4); // 3 + 2 - 1 = 4
        assert!(f.t()[0] >= 0.0);
        assert!(f.t().iter().all(|&x| x.is_finite()));
        assert!(f.y().iter().all(|&y| y.is_finite() && y >= 0.0));
      },
      other => panic!("Expected Function distribution, got {other:?}"),
    }
  }

  #[test]
  fn test_convolution_function_function_single_points() {
    let a_x = array![2.0];
    let a_y = array![3.0];
    let a = Distribution::function(a_x, a_y).unwrap();

    let b_x = array![5.0];
    let b_y = array![4.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let result = distribution_convolution(&a, &b).unwrap();
    let expected = Distribution::point(7.0, 12.0);
    assert_eq!(expected, result);
  }

  #[test]
  fn test_convolution_function_function_empty() {
    let a = Distribution::function(array![], array![]).unwrap();
    let b_x = array![1.0, 2.0];
    let b_y = array![1.0, 1.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let result = distribution_convolution(&a, &b).unwrap();
    assert_eq!(Distribution::empty(), result);
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

    let result = distribution_convolution(&a, &b).unwrap();

    match result {
      Distribution::Function(f) => {
        // Should handle different spacings correctly
        assert!(f.t().len() >= 3);
        assert!(f.t().iter().all(|&x| x.is_finite()));
        assert!(f.y().iter().all(|&y| y.is_finite() && y >= 0.0));
        // Result should span from 0+0 to 1+2 = 3
        assert!(f.t()[0] >= 0.0 - 1e-10);
        assert!(*f.t().last().unwrap() <= 3.0 + 1e-10);
      },
      other => panic!("Expected Function distribution, got {other:?}"),
    }
  }

  #[test]
  fn test_convolution_function_function_zero_width() {
    // Test point distribution convolved with function
    let a = Distribution::point(5.0, 1.0);

    let b_x = array![1.0, 2.0];
    let b_y = array![2.0, 3.0];
    let b = Distribution::function(b_x, b_y).unwrap();

    let result = distribution_convolution(&a, &b).unwrap();

    match result {
      Distribution::Function(f) => {
        // Should shift the function by the point's position
        assert!(f.t()[0] >= 6.0 - 1e-10); // 5 + 1
        assert!(*f.t().last().unwrap() <= 7.0 + 1e-10); // 5 + 2
      },
      other => panic!("Expected Function distribution, got {other:?}"),
    }
  }

  #[test]
  fn test_convolution_error_cases() {
    // Test non-uniform grids
    let non_uniform_x = array![0.0, 1.0, 3.0]; // Non-uniform spacing
    let y = array![1.0, 1.0, 1.0];
    let non_uniform = Distribution::function(non_uniform_x, y.clone()).unwrap();

    let uniform_x = array![0.0, 1.0, 2.0];
    let uniform = Distribution::function(uniform_x, y).unwrap();

    let result = distribution_convolution(&non_uniform, &uniform);
    assert!(result.is_err());
  }

  #[test]
  fn test_backward_pass_temporal_direction() {
    // Test the specific use case in backward pass: parent_time = child_time + (-branch_length)
    let child_time_dist = Distribution::point(2013.0, 1.0);
    let branch_length_dist = Distribution::point(2.5, 1.0);

    // In backward pass, we negate the branch distribution
    let negated_branch = branch_length_dist.negate();
    let parent_dist = distribution_convolution(&child_time_dist, &negated_branch).unwrap();

    match parent_dist {
      Distribution::Point(p) => {
        let parent_time = p.t();
        let expected_parent_time = 2013.0 - 2.5; // 2010.5
        assert!((parent_time - expected_parent_time).abs() < 1e-10);
        assert!(parent_time < 2013.0, "Parent should be older (earlier time) than child");
      },
      other => panic!("Expected point distribution, got {other:?}"),
    }
  }

  #[test]
  fn test_forward_pass_temporal_direction() {
    // Test the forward pass: child_time = parent_time + branch_length
    let parent_time_dist = Distribution::point(2010.0, 1.0);
    let branch_length_dist = Distribution::point(1.5, 1.0);

    let child_dist = distribution_convolution(&parent_time_dist, &branch_length_dist).unwrap();

    match child_dist {
      Distribution::Point(p) => {
        let child_time = p.t();
        let expected_child_time = 2010.0 + 1.5; // 2011.5
        assert!((child_time - expected_child_time).abs() < 1e-10);
        assert!(child_time > 2010.0, "Child should be younger (later time) than parent");
      },
      other => panic!("Expected point distribution, got {other:?}"),
    }
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

    let child_dist = distribution_convolution(&parent_dist, &branch_dist).unwrap();

    match child_dist {
      Distribution::Function(f) => {
        // Result should span from 2010+1=2011 to 2011+2=2013
        assert!(f.t()[0] >= 2011.0 - 1e-10);
        assert!(*f.t().last().unwrap() <= 2013.0 + 1e-10);
        // All probabilities should be non-negative
        assert!(f.y().iter().all(|&y| y >= 0.0));
      },
      other => panic!("Expected function distribution, got {other:?}"),
    }
  }
}
