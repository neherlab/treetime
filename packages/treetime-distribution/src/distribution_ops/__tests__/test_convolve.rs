#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::convolve::distribution_convolution;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array1, array};
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
  fn test_backward_pass_temporal_direction() -> Result<(), Report> {
    // Test the specific use case in backward pass: parent_time = child_time + (-branch_length)
    let child_time_dist = Distribution::point(2013.0, 1.0);
    let branch_length_dist = Distribution::point(2.5, 1.0);

    // In backward pass, we negate the branch distribution
    let negated_branch = branch_length_dist.negate()?;
    let actual = distribution_convolution(&child_time_dist, &negated_branch)?;

    let expected = Distribution::point(2010.5, 1.0);
    assert_eq!(expected, actual);
    Ok(())
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

  /// Regression: small dx values caused "x array must be uniformly spaced" when grid
  /// x-arrays accumulated floating-point rounding errors over many points. The fix
  /// avoids the roundtrip through explicit x-arrays by using from_start_dx_values.
  #[test]
  fn test_convolution_convolve_small_dx_function_function() -> Result<(), Report> {
    let dx = 1e-7;
    let num_points = 500;
    // Constant distributions: convolution output length = n_a + n_b - 1 = 999
    let values = Array1::from_elem(num_points, 1.0);
    let dist_a: Distribution =
      Distribution::Function(DistributionFunction::from_start_dx_values(0.0, dx, values.clone())?);
    let dist_b: Distribution = Distribution::Function(DistributionFunction::from_start_dx_values(0.0, dx, values)?);

    let actual = distribution_convolution(&dist_a, &dist_b)?;
    assert!(matches!(actual, Distribution::Function(_)));

    // Verify structural properties derived from convolution theory:
    // output length = n_a + n_b - 1, output dx preserved, output x_min = x_min_a + x_min_b
    let Distribution::Function(result_fn) = actual else {
      return Err(eyre::eyre!("expected Function variant"));
    };
    let expected_len = 2 * num_points - 1;
    assert_eq!(expected_len, result_fn.len());
    assert_abs_diff_eq!(0.0, result_fn.x_min(), epsilon = 1e-15);
    assert_abs_diff_eq!(dx, result_fn.dx(), epsilon = 1e-15);
    Ok(())
  }

  /// Regression: convolution_range_function built its result via Distribution::function(t_out, y_out)
  /// where t_out was generated from grid parameters and then validated back through Grid::from_array.
  #[test]
  fn test_convolution_convolve_small_dx_range_function() -> Result<(), Report> {
    let dx = 1e-7;
    let n = 500;
    let y = Array1::from_elem(n, 1.0);
    let func: Distribution = Distribution::Function(DistributionFunction::from_start_dx_values(0.0, dx, y)?);

    let range = Distribution::range((0.0, 2e-7), 1.0);

    let actual = distribution_convolution(&range, &func)?;
    assert!(matches!(actual, Distribution::Function(_)));

    // Range+function preserves the function's grid spacing and point count
    let Distribution::Function(f) = actual else {
      return Err(eyre::eyre!("expected Function variant"));
    };
    assert_eq!(n, f.len());
    assert_abs_diff_eq!(dx, f.dx(), epsilon = 1e-15);
    Ok(())
  }
}
