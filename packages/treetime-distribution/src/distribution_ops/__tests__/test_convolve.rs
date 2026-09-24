#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::__tests__::aliases::DistributionPlain;
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
    assert_error!(distribution_convolution(&distribution(left), &distribution(right)), expected);
  }

  #[test]
  fn test_convolution_empty() {
    let a: DistributionPlain = DistributionPlain::empty();
    let b: DistributionPlain = DistributionPlain::function(array![], array![]).unwrap();
    let actual: DistributionPlain = distribution_convolution(&a, &b).unwrap();
    let expected: DistributionPlain = DistributionPlain::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_point() {
    let a: DistributionPlain = DistributionPlain::point(2.0, 3.0);
    let b: DistributionPlain = DistributionPlain::point(5.0, 4.0);
    let actual: DistributionPlain = distribution_convolution(&a, &b).unwrap();
    let expected: DistributionPlain = DistributionPlain::point(7.0, 12.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_triangle() {
    let a = DistributionPlain::range((2.0, 4.0), 3.0);
    let b = DistributionPlain::range((6.0, 8.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![8.0, 10.0, 12.0];
      let y = array![0.0, 6.0, 0.0];
      DistributionPlain::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_trapezoid_non_uniform() {
    let a = DistributionPlain::range((2.0, 4.0), 3.0);
    let b = DistributionPlain::range((6.0, 9.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![8.0, 9.0, 10.0, 11.0, 12.0, 13.0];
      let y = array![0.0, 3.0, 6.0, 6.0, 3.0, 0.0];
      DistributionPlain::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_range_trapezoid_uniform() {
    let a = DistributionPlain::range((0.0, 2.0), 1.0);
    let b = DistributionPlain::range((3.0, 7.0), 2.0);
    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = {
      let x = array![3.0, 5.0, 7.0, 9.0];
      let y = array![0.0, 2.0, 2.0, 0.0];
      DistributionPlain::function(x, y).unwrap()
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_point_function() {
    let p = DistributionPlain::point(3.0, 2.0);

    let x = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let f = DistributionPlain::function(x, y).unwrap();
    let actual = distribution_convolution(&p, &f).unwrap();

    let x = array![3.0, 4.0, 5.0, 6.0, 7.0];
    let y = array![2.0, 4.0, 6.0, 8.0, 10.0];
    let expected = DistributionPlain::function(x, y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_function() {
    let r = DistributionPlain::range((2.0, 6.0), 0.5);

    let x = array![0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
    let y = array![0.0, 1.0, 0.0, 2.0, 1.0, 0.0];
    let f = DistributionPlain::function(x, y).unwrap();
    let actual = distribution_convolution(&r, &f).unwrap();

    let expected_t = array![4.0, 6.0, 8.0, 10.0, 12.0, 14.0];
    let expected_y = array![1.0, 1.0, 3.0, 3.0, 3.0, 1.0];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y().unwrap(), epsilon = 1e-12);
  }

  #[test]
  fn test_convolution_point_range() {
    let p = DistributionPlain::point(3.0, 2.0);
    let r = DistributionPlain::range((1.0, 4.0), 1.5);
    let actual = distribution_convolution(&p, &r).unwrap();
    let expected = DistributionPlain::range((4.0, 7.0), 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_range_point() {
    let r = DistributionPlain::range((1.0, 4.0), 1.5);
    let p = DistributionPlain::point(3.0, 2.0);
    let actual = distribution_convolution(&r, &p).unwrap();
    let expected = DistributionPlain::range((4.0, 7.0), 3.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_basic() {
    let a_x = array![0.0, 1.0, 2.0];
    let a_y = array![1.0, 2.0, 1.0];
    let a = DistributionPlain::function(a_x, a_y).unwrap();

    let b_x = array![0.0, 1.0];
    let b_y = array![1.0, 2.0];
    let b = DistributionPlain::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();

    let expected_t = array![0.0, 1.0, 2.0, 3.0];
    let expected_y = array![1.0, 4.0, 5.0, 2.0];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y().unwrap(), epsilon = 1e-12);
  }

  #[test]
  fn test_convolution_function_function_single_points() {
    let a_x = array![2.0];
    let a_y = array![3.0];
    let a = DistributionPlain::function(a_x, a_y).unwrap();

    let b_x = array![5.0];
    let b_y = array![4.0];
    let b = DistributionPlain::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = DistributionPlain::point(7.0, 12.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_function_function_empty() {
    let a = DistributionPlain::function(array![], array![]).unwrap();
    let b_x = array![1.0, 2.0];
    let b_y = array![1.0, 1.0];
    let b = DistributionPlain::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();
    let expected = DistributionPlain::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_disjoint_grids_is_not_empty() {
    let a = DistributionPlain::function(array![0.0, 1.0], array![1.0, 1.0]).unwrap();
    let b = DistributionPlain::function(array![10.0, 11.0], array![1.0, 1.0]).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();
    assert!(!matches!(actual, DistributionPlain::Empty));
  }

  #[test]
  fn test_convolution_function_function_different_spacing() {
    let a_x = array![0.0, 0.5, 1.0];
    let a_y = array![1.0, 2.0, 1.0];
    let a = DistributionPlain::function(a_x, a_y).unwrap();

    let b_x = array![0.0, 1.0, 2.0];
    let b_y = array![1.0, 1.0, 1.0];
    let b = DistributionPlain::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();

    let expected_t = array![0.0, 1.0, 2.0, 3.0];
    let expected_y = array![0.5, 2.0, 2.0, 0.5];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y().unwrap(), epsilon = 1e-12);
  }

  #[test]
  fn test_convolution_function_function_zero_width() {
    let a = DistributionPlain::point(5.0, 1.0);

    let b_x = array![1.0, 2.0];
    let b_y = array![2.0, 3.0];
    let b = DistributionPlain::function(b_x, b_y).unwrap();

    let actual = distribution_convolution(&a, &b).unwrap();

    let expected_x = array![6.0, 7.0];
    let expected_y = array![2.0, 3.0];
    let expected = DistributionPlain::function(expected_x, expected_y).unwrap();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_backward_pass_temporal_direction() -> Result<(), Report> {
    let child_time_dist = DistributionPlain::point(2013.0, 1.0);
    let branch_length_dist = DistributionPlain::point(2.5, 1.0);

    let negated_branch = branch_length_dist.negate()?;
    let actual = distribution_convolution(&child_time_dist, &negated_branch)?;

    let expected = DistributionPlain::point(2010.5, 1.0);
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_forward_pass_temporal_direction() {
    let parent_time_dist = DistributionPlain::point(2010.0, 1.0);
    let branch_length_dist = DistributionPlain::point(1.5, 1.0);

    let actual = distribution_convolution(&parent_time_dist, &branch_length_dist).unwrap();

    let expected = DistributionPlain::point(2011.5, 1.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_convolution_with_uncertainty() {
    let parent_x = array![2010.0, 2010.5, 2011.0];
    let parent_y = array![0.2, 0.6, 0.2];
    let parent_dist = DistributionPlain::function(parent_x, parent_y).unwrap();

    let branch_x = array![1.0, 1.5, 2.0];
    let branch_y = array![0.3, 0.4, 0.3];
    let branch_dist = DistributionPlain::function(branch_x, branch_y).unwrap();

    let actual = distribution_convolution(&parent_dist, &branch_dist).unwrap();

    let expected_t = array![2011.0, 2011.5, 2012.0, 2012.5, 2013.0];
    let expected_y = array![0.03, 0.13, 0.18, 0.13, 0.03];
    assert_eq!(expected_t, actual.t());
    pretty_assert_abs_diff_eq!(expected_y, actual.y().unwrap(), epsilon = 1e-12);
  }

  #[test]
  fn test_convolution_convolve_small_dx_function_function() -> Result<(), Report> {
    let dx = 1e-7;
    let num_points = 500;
    let values = Array1::from_elem(num_points, 1.0);
    let dist_a: DistributionPlain =
      DistributionPlain::Function(DistributionFunction::from_start_dx_values(0.0, dx, values.clone())?);
    let dist_b: DistributionPlain =
      DistributionPlain::Function(DistributionFunction::from_start_dx_values(0.0, dx, values)?);

    let actual = distribution_convolution(&dist_a, &dist_b)?;
    assert!(matches!(actual, DistributionPlain::Function(_)));

    let DistributionPlain::Function(result_fn) = actual else {
      return Err(eyre::eyre!("expected Function variant"));
    };
    let expected_len = 2 * num_points - 1;
    assert_eq!(expected_len, result_fn.len());
    assert_abs_diff_eq!(0.0, result_fn.x_min(), epsilon = 1e-15);
    assert_abs_diff_eq!(dx, result_fn.dx(), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_convolution_convolve_small_dx_range_function() -> Result<(), Report> {
    let dx = 1e-7;
    let n = 500;
    let y = Array1::from_elem(n, 1.0);
    let func: DistributionPlain = DistributionPlain::Function(DistributionFunction::from_start_dx_values(0.0, dx, y)?);

    let range = DistributionPlain::range((0.0, 2e-7), 1.0);

    let actual = distribution_convolution(&range, &func)?;
    assert!(matches!(actual, DistributionPlain::Function(_)));

    let DistributionPlain::Function(f) = actual else {
      return Err(eyre::eyre!("expected Function variant"));
    };
    assert_eq!(n, f.len());
    assert_abs_diff_eq!(dx, f.dx(), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_coarsen_convolution_narrow_range_keeps_fine_grid() -> Result<(), Report> {
    let fine = DistributionFunction::<f64, _>::from_start_dx_values(0.0, 0.01, array![1.0, 2.0, 1.0])?;
    assert_error!(
      Grid::from_range_dx(0.0, 0.02, 1.0),
      "Grid must have at least 2 points, got 1"
    );

    let actual: DistributionPlain = coarsen_convolution(fine.clone(), 1.0)?;

    assert_eq!(DistributionPlain::Function(fine), actual);
    Ok(())
  }

  #[test]
  fn test_coarsen_convolution_half_cell_range_resamples_to_two_points() -> Result<(), Report> {
    let fine = DistributionFunction::<f64, _>::from_start_dx_values(0.0, 0.25, array![1.0, 2.0, 1.0])?;
    let expected_points = Grid::from_range_dx(0.0, 0.5, 1.0)?.n_points();
    assert_eq!(2, expected_points);

    let DistributionPlain::Function(actual) = coarsen_convolution(fine, 1.0)? else {
      return Err(eyre::eyre!("expected Function variant"));
    };

    assert_eq!(expected_points, actual.len());
    assert_abs_diff_eq!(1.0, actual.dx(), epsilon = 1e-15);
    Ok(())
  }

  #[test]
  fn test_coarsen_convolution_wide_range_resamples_to_coarse_spacing() -> Result<(), Report> {
    let fine = DistributionFunction::<f64, _>::from_start_dx_values(0.0, 0.5, Array1::from_elem(9, 1.0))?;
    let expected_points = Grid::from_range_dx(0.0, 4.0, 1.0)?.n_points();
    assert_eq!(5, expected_points);

    let DistributionPlain::Function(actual) = coarsen_convolution(fine, 1.0)? else {
      return Err(eyre::eyre!("expected Function variant"));
    };

    assert_eq!(expected_points, actual.len());
    assert_abs_diff_eq!(1.0, actual.dx(), epsilon = 1e-15);
    Ok(())
  }

  mod helpers {
    use crate::__tests__::aliases::DistributionPlain;
    use crate::distribution_core::formula::DistributionFormula;
    use ndarray::array;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum DistributionVariant {
      Empty,
      Point,
      Range,
      Function,
      Formula,
    }

    pub(super) fn distribution(variant: DistributionVariant) -> DistributionPlain {
      match variant {
        DistributionVariant::Empty => DistributionPlain::empty(),
        DistributionVariant::Point => DistributionPlain::point(0.0, 1.0),
        DistributionVariant::Range => DistributionPlain::range((0.0, 1.0), 1.0),
        DistributionVariant::Function => {
          DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 1.0]).unwrap()
        },
        DistributionVariant::Formula => DistributionPlain::Formula(DistributionFormula::new(|_| Ok(1.0), 0.0, 1.0)),
      }
    }
  }
}
