#[cfg(test)]
mod tests {
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::policy::{NegLog, Plain};
  use crate::{Distribution, DistributionPlain};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, SoftTailLaw};
  use treetime_utils::assert_error;

  type DistFnPlain = DistributionFunction<f64, Plain>;
  type DistFnNegLog = DistributionFunction<f64, NegLog>;

  #[test]
  fn test_boundary_neglog_hard_left_returns_infinity() -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let f = f.with_left_extrap(BoundaryBehavior::Hard)?;
    // Infinity is exact -- no floating-point precision concern
    #[allow(clippy::float_cmp)]
    {
      assert_eq!(f64::INFINITY, f.interp(-1.0)?);
    }
    Ok(())
  }

  #[test]
  fn test_boundary_neglog_hard_right_returns_infinity() -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let f = f.with_right_extrap(BoundaryBehavior::Hard)?;
    #[allow(clippy::float_cmp)]
    {
      assert_eq!(f64::INFINITY, f.interp(3.0)?);
    }
    Ok(())
  }

  #[test]
  fn test_boundary_neglog_constant_allowed() -> Result<(), Report> {
    // Constant is representation-independent (it repeats a stored value), so it is allowed
    // under NegLog. The right tail returns the last stored value.
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let f = f.with_right_extrap(BoundaryBehavior::Constant)?;
    assert_ulps_eq!(2.0, f.interp(3.0)?, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_boundary_plain_hard_allowed() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let f = f.with_extrap(BoundaryBehavior::Hard)?;
    assert_ulps_eq!(0.0, f.interp(-1.0)?, max_ulps = 4);
    assert_ulps_eq!(0.0, f.interp(3.0)?, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_boundary_function_eval_errors_by_default() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let dist = Distribution::Function(f);
    assert_error!(
      dist.eval(3.0),
      "GridFn evaluated at 3.0, above the support boundary 2.0, but no extrapolation policy is set for that side"
    );
    Ok(())
  }

  #[test]
  fn test_boundary_with_extrap_noop_on_point() -> Result<(), Report> {
    // Non-Function variants have no interpolated tail, so setting a tail policy is a no-op.
    let point: DistributionPlain = Distribution::point(1.0, 2.0);
    let unchanged = point.clone().with_left_extrap(BoundaryBehavior::Hard)?;
    assert_eq!(point, unchanged);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::error_left(    (-1.0, BoundaryBehavior::Error,    BoundaryBehavior::Error),    None)]
  #[case::error_right(   ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Error),    None)]
  #[case::hard_left(     (-1.0, BoundaryBehavior::Hard,     BoundaryBehavior::Error),    None)]
  #[case::hard_right(    ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Hard),     None)]
  #[case::constant_left( (-1.0, BoundaryBehavior::Constant, BoundaryBehavior::Error),    Some(2.0))]
  #[case::constant_right(( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Constant), Some(6.0))]
  // Soft Linear tail: the function continues, so the product is non-empty. The right edge ordinate
  // is 3.0 at x_max=2.0, so the neg-log tail at t=3.0 is 3.0 + slope*(3.0-2.0) = 4.0 and the
  // product is 2.0 times it.
  #[case::linear_some_right(( 3.0, BoundaryBehavior::Error, BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 })), Some(8.0))]
  #[trace]
  fn test_boundary_multiply_point_function_outside_support(
    #[case] (t, left, right): (f64, BoundaryBehavior, BoundaryBehavior),
    #[case] expected_amplitude: Option<f64>,
  ) -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?
      .with_left_extrap(left)?
      .with_right_extrap(right)?;
    let function = Distribution::Function(f);
    let point = Distribution::point(t, 2.0);

    let actual = distribution_multiplication(&point, &function)?;
    // Oracle: kb/decisions/distribution-tails-and-arithmetic.md defines soft tails as evaluable
    // outside the grid and Hard/Error as carrying no product there. Linear is a soft tail, so it
    // extends; its value is SoftTailLaw's neg-log line y_edge + slope*(t-t_edge).
    let expected = expected_amplitude.map_or_else(Distribution::empty, |amplitude| Distribution::point(t, amplitude));
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_boundary_multiply_composes_linear_tail_slopes() -> Result<(), Report> {
    // Two functions on the same support, each with a right soft tail. Multiplication is addition
    // in neg-log space, so the product's right tail decays with the sum of the two slopes,
    // anchored on the product's own right edge value.
    let (slope_a, slope_b) = (0.6, 0.9);
    let fa: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![4.0, 3.0, 2.0])?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: slope_a }))?;
    let fb: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 1.5, 2.5])?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: slope_b }))?;

    let product = distribution_multiplication(&Distribution::Function(fa), &Distribution::Function(fb))?;

    // The product's right edge ordinate is 2.0 * 2.5 = 5.0 at x_max = 2.0. Evaluated 0.5 beyond the
    // edge, the composed neg-log tail is 5.0 + (slope_a + slope_b) * 0.5.
    let expected = 5.0 + (slope_a + slope_b) * 0.5;
    assert_ulps_eq!(expected, product.eval(2.5)?, max_ulps = 8);
    Ok(())
  }

  #[test]
  fn test_boundary_multiply_point_inside_function_is_point() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let function = Distribution::Function(f);
    let point_inside = Distribution::point(1.0, 2.0);
    let actual = distribution_multiplication(&point_inside, &function)?;
    // function(1.0) = 2.0, point amplitude 2.0, product 4.0 at t = 1.0
    assert_eq!(Distribution::point(1.0, 4.0), actual);
    Ok(())
  }
}
