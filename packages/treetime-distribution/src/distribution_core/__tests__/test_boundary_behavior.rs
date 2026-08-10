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
  use treetime_grid::BoundaryBehavior;
  use treetime_utils::assert_error;

  type DistFnPlain = DistributionFunction<f64, Plain>;
  type DistFnNegLog = DistributionFunction<f64, NegLog>;

  #[test]
  fn test_boundary_neglog_hard_left_rejected() -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    assert_error!(
      f.with_left_extrap(BoundaryBehavior::Hard(None)),
      "Refusing a Hard boundary tail: it writes 0.0 outside support, which is the multiplicative identity (probability one), not zero probability, under this distribution's negative-log representation"
    );
    Ok(())
  }

  #[test]
  fn test_boundary_neglog_hard_right_rejected() -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    assert_error!(
      f.with_right_extrap(BoundaryBehavior::Hard(None)),
      "Refusing a Hard boundary tail: it writes 0.0 outside support, which is the multiplicative identity (probability one), not zero probability, under this distribution's negative-log representation"
    );
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
    let f = f.with_extrap(BoundaryBehavior::Hard(None))?;
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
    // Non-Function variants have no interpolated tail, so setting a tail policy is a no-op
    // and never rejects (even the Hard tail that a Function would reject under NegLog).
    let point: DistributionPlain = Distribution::point(1.0, 2.0);
    let unchanged = point.clone().with_left_extrap(BoundaryBehavior::Hard(None))?;
    assert_eq!(point, unchanged);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::error_left(    (-1.0, BoundaryBehavior::Error,    BoundaryBehavior::Error),    None)]
  #[case::error_right(   ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Error),    None)]
  #[case::hard_left(     (-1.0, BoundaryBehavior::Hard(None),     BoundaryBehavior::Error),    None)]
  #[case::hard_right(    ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Hard(None)),     None)]
  #[case::constant_left( (-1.0, BoundaryBehavior::Constant, BoundaryBehavior::Error),    Some(2.0))]
  #[case::constant_right(( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Constant), Some(6.0))]
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
    // Oracle: kb/decisions/distribution-tails-and-arithmetic.md#L74-L84 defines Constant
    // as evaluable outside the grid and Hard/Error as carrying no product there.
    let expected = expected_amplitude.map_or_else(Distribution::empty, |amplitude| Distribution::point(t, amplitude));
    assert_eq!(expected, actual);
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
