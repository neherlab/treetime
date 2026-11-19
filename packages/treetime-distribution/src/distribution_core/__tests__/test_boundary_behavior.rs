#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionPlain;
  use crate::Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use eyre::Report;
  use ndarray::array;
  use rstest::rstest;
  use treetime_grid::MaxGridPoints;
  use treetime_grid::{BoundaryBehavior, Side, SoftTailLaw};
  use treetime_grid::{NegLog, Plain};
  use treetime_utils::assert_error;
  use treetime_utils::pretty_assert_ulps_eq;

  type DistFnPlain = DistributionFunction<f64, Plain>;
  type DistFnNegLog = DistributionFunction<f64, NegLog>;

  #[rustfmt::skip]
  #[rstest]
  #[case::left( Side::Left,  -1.0)]
  #[case::right(Side::Right,  3.0)]
  #[trace]
  fn test_boundary_neglog_hard_side_is_zero_probability_outside_support(
    #[case] side: Side,
    #[case] t: f64,
  ) -> Result<(), Report> {
    let f: DistFnNegLog = DistributionFunction::from_range_values((0.0, 2.0), array![0.0, 1.0, 2.0])?;
    let f = match side {
      Side::Left => f.with_left_extrap(BoundaryBehavior::Hard)?,
      Side::Right => f.with_right_extrap(BoundaryBehavior::Hard)?,
    };
    assert_eq!(f64::INFINITY, f.interp(t)?);
    Ok(())
  }

  #[test]
  fn test_boundary_plain_hard_sides_are_zero_outside_support() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let f = f
      .with_left_extrap(BoundaryBehavior::Hard)?
      .with_right_extrap(BoundaryBehavior::Hard)?;
    pretty_assert_ulps_eq!(0.0, f.interp(-1.0)?, max_ulps = 4);
    pretty_assert_ulps_eq!(0.0, f.interp(3.0)?, max_ulps = 4);
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

    pretty_assert_ulps_eq!(0.0, empty.eval(1.0)?, max_ulps = 4);
    pretty_assert_ulps_eq!(0.0, point.eval(0.0)?, max_ulps = 4);
    assert!(range.eval(3.0)?.is_infinite());
    Ok(())
  }

  #[test]
  fn test_boundary_formula_ignores_soft_tail_fitting() -> Result<(), Report> {
    let formula: Distribution<NegLog> = Distribution::Formula(DistributionFormula::new(Ok, 0.0, 2.0));
    let fitted = formula.fit_soft_tail(Side::Right, 10)?;

    assert!(matches!(fitted, Distribution::Formula(_)));
    pretty_assert_ulps_eq!(2.5, fitted.eval(2.5)?, max_ulps = 8);
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
    pretty_assert_ulps_eq!(2.5, distribution.eval(2.5)?, max_ulps = 8);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::error_left(    (-1.0, BoundaryBehavior::Error,    BoundaryBehavior::Error))]
  #[case::error_right(   ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Error))]
  #[case::hard_left(     (-1.0, BoundaryBehavior::Hard,     BoundaryBehavior::Error))]
  #[case::hard_right(    ( 3.0, BoundaryBehavior::Error,    BoundaryBehavior::Hard))]
  #[trace]
  fn test_boundary_multiply_point_function_outside_hard_support_is_empty(
    #[case] (t, left, right): (f64, BoundaryBehavior, BoundaryBehavior),
  ) -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?
      .with_left_extrap(left)?
      .with_right_extrap(right)?;
    let function = Distribution::Function(f);
    let point = Distribution::point(t, 2.0);

    let actual = distribution_multiplication(&point, &function, MaxGridPoints::default())?;
    assert_eq!(Distribution::empty(), actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::left( (-1.0, BoundaryBehavior::Linear(SoftTailLaw { slope: -1.0 }), BoundaryBehavior::Error), 2.0 * (-1.0_f64).exp())]
  #[case::right(( 3.0, BoundaryBehavior::Error, BoundaryBehavior::Linear(SoftTailLaw { slope: 1.0 })), 2.0 * 3.0 * (-1.0_f64).exp())]
  #[trace]
  fn test_boundary_multiply_point_function_in_soft_tail_reads_the_tail(
    #[case] (t, left, right): (f64, BoundaryBehavior, BoundaryBehavior),
    #[case] expected_amplitude: f64,
  ) -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?
      .with_left_extrap(left)?
      .with_right_extrap(right)?;
    let function = Distribution::Function(f);
    let point = Distribution::point(t, 2.0);

    let actual = distribution_multiplication(&point, &function, MaxGridPoints::default())?;
    let DistributionPlain::Point(actual) = actual else {
      panic!("expected a point, got {actual:?}");
    };
    pretty_assert_ulps_eq!(t, actual.t(), max_ulps = 0);
    pretty_assert_ulps_eq!(
      expected_amplitude,
      actual.amplitude(),
      max_ulps = 4,
      "point amplitude 2 times the edge density e^0 or 3, decayed by e^-1 over one unit at slope 1 in neg-log. Bound: one ulp for each of ln 3, the subtraction and exp in the evaluation, and exp in the oracle"
    );
    Ok(())
  }

  #[test]
  fn test_boundary_multiply_composes_linear_tail_slopes() -> Result<(), Report> {
    let (slope_a, slope_b) = (0.6, 0.9);
    let fa: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![4.0, 3.0, 2.0])?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: slope_a }))?;
    let fb: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 1.5, 2.5])?
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: slope_b }))?;

    let product = distribution_multiplication(
      &Distribution::Function(fa),
      &Distribution::Function(fb),
      MaxGridPoints::default(),
    )?;

    let expected = 2.0 * (-slope_a * 0.5).exp() * 2.5 * (-slope_b * 0.5).exp();
    pretty_assert_ulps_eq!(
      expected,
      product.eval(2.5)?,
      max_ulps = 11,
      "product of the operand tails at 2.5: 2 e^(-0.6 * 0.5) times 2.5 e^(-0.9 * 0.5). Bound: one ulp for each of the 5 roundings and libm calls of the evaluation (product, ln, slope times step, sum, exp) and the 6 of the oracle"
    );
    Ok(())
  }

  #[test]
  fn test_boundary_multiply_point_inside_function_is_point() -> Result<(), Report> {
    let f: DistFnPlain = DistributionFunction::from_range_values((0.0, 2.0), array![1.0, 2.0, 3.0])?;
    let function = Distribution::Function(f);
    let point_inside = Distribution::point(1.0, 2.0);
    let actual = distribution_multiplication(&point_inside, &function, MaxGridPoints::default())?;
    assert_eq!(Distribution::point(1.0, 4.0), actual);
    Ok(())
  }
}
