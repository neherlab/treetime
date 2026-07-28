#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_ops::divide::distribution_division;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::BoundaryBehavior;
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
    // Oracle: kb/decisions/timetree-inference-pass-boundary-tails.md, "Scope and interaction".
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
