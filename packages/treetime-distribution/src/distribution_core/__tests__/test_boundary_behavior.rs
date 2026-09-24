#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionPlain;
  use crate::Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::policy::{NegLog, Plain};
  use approx::assert_ulps_eq;
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, Side, SoftTailLaw};
  use treetime_utils::assert_error;

  type DistFnPlain = DistributionFunction<f64, Plain>;
  type DistFnNegLog = DistributionFunction<f64, NegLog>;

  #[test]
  fn test_boundary_neglog_hard_left_returns_infinity() -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let f = f.with_left_extrap(BoundaryBehavior::Hard)?;
    #[expect(clippy::float_cmp, reason = "hard boundaries return exactly infinity")]
    {
      assert_eq!(f64::INFINITY, f.interp(-1.0)?);
    }
    Ok(())
  }

  #[test]
  fn test_boundary_neglog_hard_right_returns_infinity() -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let f = f.with_right_extrap(BoundaryBehavior::Hard)?;
    #[expect(clippy::float_cmp, reason = "hard boundaries return exactly infinity")]
    {
      assert_eq!(f64::INFINITY, f.interp(3.0)?);
    }
    Ok(())
  }

  #[test]
  fn test_boundary_plain_hard_allowed() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let f = f
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
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
  fn test_boundary_compact_variant_remains_hard() -> Result<(), Report> {
    let point: DistributionPlain = Distribution::point(1.0, 2.0);
    let unchanged = point
      .clone()
      .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }))?;
    assert_eq!(point, unchanged);
    Ok(())
  }

  #[test]
  fn test_boundary_compact_variants_evaluate_to_probability_zero_outside_support() -> Result<(), Report> {
    let empty: DistributionPlain = Distribution::empty();
    let point: DistributionPlain = Distribution::point(1.0, 2.0);
    let range: Distribution<NegLog> = Distribution::range((0.0, 2.0), 3.0);

    assert_ulps_eq!(0.0, empty.eval(1.0)?, max_ulps = 4);
    assert_ulps_eq!(0.0, point.eval(0.0)?, max_ulps = 4);
    assert!(range.eval(3.0)?.is_infinite());
    Ok(())
  }

  #[test]
  fn test_boundary_formula_ignores_soft_tail_fitting() -> Result<(), Report> {
    let formula: Distribution<NegLog> = Distribution::Formula(DistributionFormula::new(Ok, 0.0, 2.0));
    let fitted = formula.fit_soft_tail(Side::Right, 10)?;

    assert!(matches!(fitted, Distribution::Formula(_)));
    assert_ulps_eq!(2.5, fitted.eval(2.5)?, max_ulps = 8);
    Ok(())
  }

  #[test]
  fn test_boundary_distribution_function_soft_tail_fitting() -> Result<(), Report> {
    let function: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let distribution = Distribution::Function(function).fit_soft_tail(Side::Right, 3)?;

    let Distribution::Function(fitted) = &distribution else {
      panic!("expected a function distribution, got {distribution:?}");
    };
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 }),
      fitted.right_extrap()
    );
    assert_ulps_eq!(2.5, distribution.eval(2.5)?, max_ulps = 8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::error_left(    (-1.0, BoundaryBehavior::Error,    BoundaryBehavior::Error),    None)]
  #[case::error_right(   ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Error),    None)]
  #[case::hard_left(     (-1.0, BoundaryBehavior::Hard,     BoundaryBehavior::Error),    None)]
  #[case::hard_right(    ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Hard),     None)]
  #[case::linear_some_left( (-1.0, BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }), BoundaryBehavior::Error), Some(4.0))]
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
    let expected = expected_amplitude.map_or_else(Distribution::empty, |amplitude| Distribution::point(t, amplitude));
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_boundary_multiply_composes_linear_tail_slopes() -> Result<(), Report> {
    let (slope_a, slope_b) = (0.6, 0.9);
    let fa: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![4.0, 3.0, 2.0])?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: slope_a }))?;
    let fb: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 1.5, 2.5])?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: slope_b }))?;

    let product = distribution_multiplication(&Distribution::Function(fa), &Distribution::Function(fb))?;

    let expected = 5.0 + f64::midpoint(slope_a, slope_b);
    assert_ulps_eq!(expected, product.eval(2.5)?, max_ulps = 8);
    Ok(())
  }

  #[test]
  fn test_boundary_multiply_point_inside_function_is_point() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let function = Distribution::Function(f);
    let point_inside = Distribution::point(1.0, 2.0);
    let actual = distribution_multiplication(&point_inside, &function)?;
    assert_eq!(Distribution::point(1.0, 4.0), actual);
    Ok(())
  }
}
