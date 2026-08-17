#[cfg(test)]
mod tests {
  use crate::DistributionNegLog;
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::divide::distribution_division;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  use self::helpers::{DistributionVariant, distribution};

  const TINY_NUMBER: f64 = 1e-10;

  #[rustfmt::skip]
  #[rstest]
  #[case::formula_empty(   DistributionVariant::Formula,  DistributionVariant::Empty,    "Cannot divide Formula by Empty: operation not implemented")]
  #[case::formula_point(   DistributionVariant::Formula,  DistributionVariant::Point,    "Cannot divide Formula by Point: operation not implemented")]
  #[case::formula_range(   DistributionVariant::Formula,  DistributionVariant::Range,    "Cannot divide Formula by Range: operation not implemented")]
  #[case::formula_function(DistributionVariant::Formula,  DistributionVariant::Function, "Cannot divide Formula by Function: operation not implemented")]
  #[case::formula_formula( DistributionVariant::Formula,  DistributionVariant::Formula,  "Cannot divide Formula by Formula: operation not implemented")]
  #[case::empty_formula(   DistributionVariant::Empty,    DistributionVariant::Formula,  "Cannot divide Empty by Formula: operation not implemented")]
  #[case::point_formula(   DistributionVariant::Point,    DistributionVariant::Formula,  "Cannot divide Point by Formula: operation not implemented")]
  #[case::range_formula(   DistributionVariant::Range,    DistributionVariant::Formula,  "Cannot divide Range by Formula: operation not implemented")]
  #[case::function_formula(DistributionVariant::Function, DistributionVariant::Formula,  "Cannot divide Function by Formula: operation not implemented")]
  #[trace]
  fn test_divide_formula_combinations_return_errors(
    #[case] dividend: DistributionVariant,
    #[case] divisor: DistributionVariant,
    #[case] expected: &str,
  ) {
    // Oracle: kb/issues/H-distribution-result-api-panics-on-formula.md.
    assert_error!(
      distribution_division(&distribution(dividend), &distribution(divisor)),
      expected
    );
  }

  #[test]
  fn test_divide_empty_by_any() {
    let empty = Distribution::empty();
    let point = Distribution::point(1.0, 2.0);
    let actual = distribution_division(&empty, &point).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_by_empty_fails() {
    let point = Distribution::point(1.0, 2.0);
    let empty = Distribution::empty();
    assert_error!(
      distribution_division(&point, &empty),
      "Cannot divide by empty distribution"
    );
  }

  #[test]
  fn test_divide_point_by_function() {
    let point = Distribution::point(2.0, 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&point, &func).unwrap();
    let expected = Distribution::point(2.0, 2.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected_y = array![5.0, 5.0, 6.0, 5.0, 5.0];
    let expected = Distribution::function(t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_by_zero_handled() {
    let t = array![0.0, 1.0, 2.0];
    let y1 = array![10.0, 20.0, 30.0];
    let y2 = array![2.0, 0.0, 5.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected_y = array![5.0, 20.0 / TINY_NUMBER, 6.0];
    let expected = Distribution::function(t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_full_overlap() {
    // Range [1.0, 3.0] overlaps with function at points [1.0, 2.0, 3.0]
    let range = Distribution::range((1.0, 3.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();
    // Result should only cover the overlap region [1.0, 3.0]
    let expected = Distribution::function(array![1.0, 2.0, 3.0], array![5.0, 2.0, 2.5]).unwrap();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::contained(        (1.0, 3.0), vec![1.0, 2.0, 3.0],                         vec![6.0, 4.0, 3.0])]
  #[case::left_partial(    (-1.0, 2.5), vec![0.0, 5.0 / 6.0, 5.0 / 3.0, 2.5],       vec![12.0, 72.0 / 11.0, 4.5, 24.0 / 7.0])]
  #[case::right_partial(    (1.5, 5.0), vec![1.5, 7.0 / 3.0, 19.0 / 6.0, 4.0],      vec![4.8, 3.6, 72.0 / 25.0, 2.4])]
  #[case::no_interior_knot( (1.2, 1.8), vec![1.2, 1.8],                              vec![60.0 / 11.0, 30.0 / 7.0])]
  #[trace]
  fn test_divide_range_by_function_preserves_analytical_overlap(
    #[case] range_bounds: (f64, f64),
    #[case] expected_t: Vec<f64>,
    #[case] expected_y: Vec<f64>,
  ) {
    let range = Distribution::range(range_bounds, 12.0);
    let function = Distribution::function(
      array![0.0, 1.0, 2.0, 3.0, 4.0],
      array![1.0, 2.0, 3.0, 4.0, 5.0],
    )
    .unwrap();

    let actual = distribution_division(&range, &function).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };
    pretty_assert_ulps_eq!(Array1::from_vec(expected_t), actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(Array1::from_vec(expected_y), actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_divide_range_by_function_endpoint_contact_returns_point() {
    let range = Distribution::range((2.0, 3.0), 12.0);
    let function = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_division(&range, &function).unwrap();
    // Oracle: v0 `Distribution.divide()` converts one surviving endpoint knot to a delta.
    let expected = Distribution::point(2.0, 4.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_no_overlap() {
    let range = Distribution::range((5.0, 6.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0];
    let y = array![1.0, 2.0, 4.0, 8.0];
    let func = Distribution::function(t, y).unwrap();

    let actual = distribution_division(&range, &func).unwrap();
    let expected = Distribution::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_same_grid() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = Distribution::function(t.clone(), y1).unwrap();
    let divisor = Distribution::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected = Distribution::function(t, array![5.0, 5.0, 6.0, 5.0, 5.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_different_grids() {
    let t1 = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let dividend = Distribution::function(t1.clone(), y1).unwrap();

    let t2 = array![0.0, 2.0, 4.0];
    let y2 = array![2.0, 5.0, 10.0];
    let divisor = Distribution::function(t2, y2).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let expected = Distribution::function(t1, array![5.0, 20.0 / 3.5, 6.0, 40.0 / 7.5, 5.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_neglog_bounds_quotient_to_hard_divisor() {
    // A NegLog dividend wider than its hard-bounded divisor. Under NegLog a hard boundary reads +inf
    // (zero probability) beyond the edge, so extending the quotient into the divisor's tail would
    // compute `dividend - (+inf) = -inf` and spike the result, collapsing the downstream forward
    // message. The quotient must be bounded to the divisor's real grid support and stay finite.
    // Oracle: kb/decisions/distribution-tails-and-arithmetic.md (Division).
    let dividend =
      DistributionNegLog::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![0.4, 0.2, 0.0, 0.2, 0.4]).unwrap();
    let divisor = DistributionNegLog::function(array![1.0, 2.0, 3.0], array![0.1, 0.0, 0.1])
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();

    let crate::Distribution::Function(function) = actual else {
      panic!("expected Function, got {actual:?}");
    };
    pretty_assert_ulps_eq!(1.0, function.x_min(), max_ulps = 4);
    pretty_assert_ulps_eq!(3.0, function.x_max(), max_ulps = 4);
    assert!(
      function.y().iter().all(|value| value.is_finite()),
      "quotient must stay finite, got {:?}",
      function.y()
    );
  }

  #[test]
  fn test_divide_function_by_function_uses_exact_intersection() {
    let dividend = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 4.0, 6.0, 8.0, 10.0]).unwrap();
    let divisor = Distribution::function(array![1.5, 2.5, 3.5], array![2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();
    let expected = Distribution::function(array![1.5, 2.5, 3.5], array![2.5, 3.5, 4.5]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_honors_explicit_divisor_tails() {
    let dividend =
      Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![10.0, 20.0, 30.0, 40.0, 50.0]).unwrap();
    let divisor = Distribution::function(array![1.0, 2.0, 3.0], array![2.0, 2.0, 2.0])
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Constant)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Constant)
      .unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();
    // Oracle: kb/decisions/distribution-tails-and-arithmetic.md, division tail rules.
    let expected =
      Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![5.0, 10.0, 15.0, 20.0, 25.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_endpoint_contact_returns_point() {
    let dividend = Distribution::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let divisor = Distribution::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_division(&dividend, &divisor).unwrap();
    // Oracle: v0 `Distribution.divide()` converts one surviving endpoint knot to a delta.
    let expected = Distribution::point(1.0, 0.6);
    assert_eq!(expected, actual);
  }

  /// The quotient inherits the dividend's per-side tail when the dividend binds the intersection.
  ///
  /// The dividend spans `[2, 10]` with a soft `Linear` left law and a `Hard` right edge; the divisor
  /// spans the wider `[0, 12]` with default `Error` tails, so its edges never truncate inside the
  /// dividend and the intersection is exactly `[2, 10]`. The quotient therefore carries the
  /// dividend's Linear-left / Hard-right policy, not the `Error` a bare rebuild would leave.
  #[test]
  fn test_divide_function_by_function_inherits_dividend_tails() {
    let left = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.5 });
    let dividend =
      DistributionFunction::from_start_dx_values(2.0, 1.0, array![9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        .unwrap()
        .with_left_extrap(left)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap();
    let divisor = DistributionFunction::from_start_dx_values(0.0, 1.0, Array1::from_elem(13, 2.0)).unwrap();

    let actual = distribution_division(&Distribution::Function(dividend), &Distribution::Function(divisor)).unwrap();

    let Distribution::Function(f) = actual else {
      panic!("expected a Function quotient");
    };
    assert_eq!(left, f.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
  }

  /// A divisor `Error` bound strictly inside the dividend truncates the result grid, so that side is
  /// `Error` (the quotient is undefined past the divisor's own edge).
  ///
  /// The dividend spans `[2, 10]` (Linear left, Hard right); the divisor spans `[4, 8]` with `Error`
  /// tails, so it truncates both sides inward of the dividend and the intersection is `[4, 8]`. Both
  /// result sides are `Error`, overriding the dividend's inherited tails.
  #[test]
  fn test_divide_function_by_function_divisor_error_truncation_yields_error_tails() {
    let dividend =
      DistributionFunction::from_start_dx_values(2.0, 1.0, array![9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        .unwrap()
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.5 }))
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap();
    let divisor = DistributionFunction::from_start_dx_values(4.0, 1.0, array![2.0, 2.0, 2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_division(&Distribution::Function(dividend), &Distribution::Function(divisor)).unwrap();

    let Distribution::Function(f) = actual else {
      panic!("expected a Function quotient");
    };
    assert_eq!(BoundaryBehavior::Error, f.left_extrap());
    assert_eq!(BoundaryBehavior::Error, f.right_extrap());
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
