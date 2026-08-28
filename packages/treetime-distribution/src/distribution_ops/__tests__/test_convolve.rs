#[cfg(test)]
mod tests {
  use crate::DistributionNegLog;
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::convolve::{coarsen_convolution, distribution_convolution};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use ndarray::{Array1, array};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_grid::grid::Grid;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  use self::helpers::{DistributionVariant, distribution};

  /// Analytic oracle: convolving two unit-variance Gaussians yields a Gaussian with variance
  /// `sigma1^2 + sigma2^2 = 2` and mean `mu1 + mu2 = 0`. In NegLog the ordinate is the negative
  /// log likelihood `y(t) = (t - mu)^2 / (2 sigma^2)`, so the result must peak at `t = 0` with
  /// near-peak curvature `1 / (2 * 2)`: `y(1) - y(0) = 0.25` and `y(2) - y(0) = 1.0`.
  ///
  /// Exercises the plain-space FFT round-trip in the trusted bulk. The far tails (reconstructed by
  /// log-linear fit) need separate validation once this path is wired into the pipeline.
  #[test]
  fn test_convolve_neglog_gaussian_variances_add() -> Result<(), Report> {
    let gaussian = |sigma: f64| -> Result<DistributionNegLog, Report> {
      let t = Array1::linspace(-6.0, 6.0, 241);
      let y = t.mapv(|t| t * t / (2.0 * sigma * sigma));
      DistributionNegLog::function(t, y)
    };

    let result = distribution_convolution(&gaussian(1.0)?, &gaussian(1.0)?)?;

    let y0 = result.eval(0.0)?;
    assert_abs_diff_eq!(0.0, result.likely_time().unwrap(), epsilon = 0.05);
    assert_abs_diff_eq!(0.25, result.eval(1.0)? - y0, epsilon = 1e-3);
    assert_abs_diff_eq!(1.00, result.eval(2.0)? - y0, epsilon = 1e-3);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::formula_empty(   DistributionVariant::Formula,  DistributionVariant::Empty,    "Cannot convolve formula with empty: operation not implemented")]
  #[case::formula_point(   DistributionVariant::Formula,  DistributionVariant::Point,    "Cannot convolve formula with point: operation not implemented")]
  #[case::formula_range(   DistributionVariant::Formula,  DistributionVariant::Range,    "Cannot convolve formula with range: operation not implemented")]
  #[case::formula_function(DistributionVariant::Formula,  DistributionVariant::Function, "Cannot convolve formula with function: operation not implemented")]
  #[case::formula_formula( DistributionVariant::Formula,  DistributionVariant::Formula,  "Cannot convolve formula with formula: operation not implemented")]
  #[case::empty_formula(   DistributionVariant::Empty,    DistributionVariant::Formula,  "Cannot convolve empty with formula: operation not implemented")]
  #[case::point_formula(   DistributionVariant::Point,    DistributionVariant::Formula,  "Cannot convolve point with formula: operation not implemented")]
  #[case::range_formula(   DistributionVariant::Range,    DistributionVariant::Formula,  "Cannot convolve range with formula: operation not implemented")]
  #[case::function_formula(DistributionVariant::Function, DistributionVariant::Formula,  "Cannot convolve function with formula: operation not implemented")]
  #[trace]
  fn test_convolve_formula_combinations_return_errors(
    #[case] left: DistributionVariant,
    #[case] right: DistributionVariant,
    #[case] expected: &str,
  ) {
    // Oracle: kb/issues/H-distribution-result-api-panics-on-formula.md.
    assert_error!(distribution_convolution(&distribution(left), &distribution(right)), expected);
  }

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

    // Convolution routes through a shared log-space peak-normalization, so the plain result carries
    // ULP-level round-trip noise. Compare the grid exactly and the amplitudes at a tight tolerance.
    let expected_t = array![4.0, 6.0, 8.0, 10.0, 12.0, 14.0];
    let expected_y = array![1.0, 1.0, 3.0, 3.0, 3.0, 1.0];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y(), epsilon = 1e-12);
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

    // Shared log-space peak-normalization adds ULP-level round-trip noise to the plain result.
    let expected_t = array![0.0, 1.0, 2.0, 3.0];
    let expected_y = array![1.0, 4.0, 5.0, 2.0];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y(), epsilon = 1e-12);
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

  /// Convolution is a Minkowski sum, so two mass-bearing operands on disjoint grids still produce a
  /// non-empty result whose support is the sum of the grids. This is why the empty-result guard
  /// models a mass-bearing convolution operand as unbounded: only an actually empty operand may
  /// yield `Empty`, never disjoint grids (which would be legitimate empties under multiplication).
  #[test]
  fn test_convolution_disjoint_grids_is_not_empty() {
    let a = Distribution::function(array![0.0, 1.0], array![1.0, 1.0]).unwrap();
    let b = Distribution::function(array![10.0, 11.0], array![1.0, 1.0]).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();
    assert!(!matches!(actual, Distribution::Empty));
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

    // Shared log-space peak-normalization and the plain-space FFT add ULP-level round-trip noise to
    // the plain result. Compare the grid exactly and the amplitudes at a tight tolerance.
    let expected_t = array![0.0, 1.0, 2.0, 3.0];
    let expected_y = array![0.5, 2.0, 2.0, 0.5];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y(), epsilon = 1e-12);
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

    // Shared log-space peak-normalization adds ULP-level round-trip noise to the plain result.
    let expected_t = array![2011.0, 2011.5, 2012.0, 2012.5, 2013.0];
    let expected_y = array![0.03, 0.13, 0.18, 0.13, 0.03];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y(), epsilon = 1e-12);
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

  /// The fine-grid convolution result is coarsened to the coarser operand spacing. When the trusted
  /// bulk is narrower than half a coarse cell, `Grid::from_range_dx` cannot build a `>= 2`-point
  /// coarse grid, so `coarsen_convolution` must keep the fine-grid result unchanged rather than fail.
  ///
  /// Oracle: `Grid::from_range_dx(0.0, 0.02, 1.0)` rejects a single-point grid, so the coarse
  /// spacing is unrepresentable and the fine result is the only valid output.
  #[test]
  fn test_coarsen_convolution_narrow_range_keeps_fine_grid() -> Result<(), Report> {
    let fine = DistributionFunction::<f64, _>::from_start_dx_values(0.0, 0.01, array![1.0, 2.0, 1.0])?;
    assert_error!(
      Grid::from_range_dx(0.0, 0.02, 1.0),
      "Grid must have at least 2 points, got 1"
    );

    let actual: Distribution = coarsen_convolution(fine.clone(), 1.0)?;

    assert_eq!(Distribution::Function(fine), actual);
    Ok(())
  }

  /// A trusted bulk exactly `0.5 * coarse_dx` wide rounds up to a two-point coarse grid, so it
  /// coarsens rather than staying fine. Oracle: `Grid::from_range_dx(0.0, 0.5, 1.0)` yields two points.
  #[test]
  fn test_coarsen_convolution_half_cell_range_resamples_to_two_points() -> Result<(), Report> {
    let fine = DistributionFunction::<f64, _>::from_start_dx_values(0.0, 0.25, array![1.0, 2.0, 1.0])?;
    let expected_points = Grid::from_range_dx(0.0, 0.5, 1.0)?.n_points();
    assert_eq!(2, expected_points);

    let Distribution::Function(actual) = coarsen_convolution(fine, 1.0)? else {
      return Err(eyre::eyre!("expected Function variant"));
    };

    assert_eq!(expected_points, actual.len());
    assert_abs_diff_eq!(1.0, actual.dx(), epsilon = 1e-15);
    Ok(())
  }

  /// A wide trusted bulk coarsens to the coarse operand spacing. Oracle: `Grid::from_range_dx(0.0,
  /// 4.0, 1.0)` yields five points at spacing `1.0`, the grid `resample_dx(1.0)` targets.
  #[test]
  fn test_coarsen_convolution_wide_range_resamples_to_coarse_spacing() -> Result<(), Report> {
    let fine = DistributionFunction::<f64, _>::from_start_dx_values(0.0, 0.5, Array1::from_elem(9, 1.0))?;
    let expected_points = Grid::from_range_dx(0.0, 4.0, 1.0)?.n_points();
    assert_eq!(5, expected_points);

    let Distribution::Function(actual) = coarsen_convolution(fine, 1.0)? else {
      return Err(eyre::eyre!("expected Function variant"));
    };

    assert_eq!(expected_points, actual.len());
    assert_abs_diff_eq!(1.0, actual.dx(), epsilon = 1e-15);
    Ok(())
  }

  mod helpers {
    use crate::DistributionPlain as Distribution;
    use crate::distribution_core::formula::DistributionFormula;
    use ndarray::array;

    #[derive(Clone, Copy, Debug)]
    pub enum DistributionVariant {
      Empty,
      Point,
      Range,
      Function,
      Formula,
    }

    pub fn distribution(variant: DistributionVariant) -> Distribution {
      match variant {
        DistributionVariant::Empty => Distribution::empty(),
        DistributionVariant::Point => Distribution::point(0.0, 1.0),
        DistributionVariant::Range => Distribution::range((0.0, 1.0), 1.0),
        DistributionVariant::Function => Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]).unwrap(),
        DistributionVariant::Formula => Distribution::Formula(DistributionFormula::new(|_| Ok(1.0), 0.0, 1.0)),
      }
    }
  }
}
