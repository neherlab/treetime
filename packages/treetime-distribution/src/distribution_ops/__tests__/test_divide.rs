#[cfg(test)]
mod tests {
  use self::helpers::{
    DistributionVariant, distribution, end_index, end_tail, hard_end, power_law_to_hard_end, power_law_values,
  };
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::__tests__::aliases::DistributionPlain;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::divide::distribution_division;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::policy::NegLog;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::MaxGridPoints;
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, HardApproachLaw, Side, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  const TINY_NUMBER: f64 = 1e-10;

  #[rustfmt::skip]
  #[rstest]
  #[case::formula_empty(   DistributionVariant::Formula,  DistributionVariant::Empty,    "Cannot divide formula by empty: operation not implemented")]
  #[case::formula_point(   DistributionVariant::Formula,  DistributionVariant::Point,    "Cannot divide formula by point: operation not implemented")]
  #[case::formula_range(   DistributionVariant::Formula,  DistributionVariant::Range,    "Cannot divide formula by range: operation not implemented")]
  #[case::formula_function(DistributionVariant::Formula,  DistributionVariant::Function, "Cannot divide formula by function: operation not implemented")]
  #[case::formula_formula( DistributionVariant::Formula,  DistributionVariant::Formula,  "Cannot divide formula by formula: operation not implemented")]
  #[case::empty_formula(   DistributionVariant::Empty,    DistributionVariant::Formula,  "Cannot divide empty by formula: operation not implemented")]
  #[case::point_formula(   DistributionVariant::Point,    DistributionVariant::Formula,  "Cannot divide point by formula: operation not implemented")]
  #[case::range_formula(   DistributionVariant::Range,    DistributionVariant::Formula,  "Cannot divide range by formula: operation not implemented")]
  #[case::function_formula(DistributionVariant::Function, DistributionVariant::Formula,  "Cannot divide function by formula: operation not implemented")]
  #[trace]
  fn test_divide_formula_combinations_return_errors(
    #[case] dividend: DistributionVariant,
    #[case] divisor: DistributionVariant,
    #[case] expected: &str,
  ) {
    assert_error!(
      distribution_division(&distribution(dividend), &distribution(divisor), MaxGridPoints::default()),
      expected
    );
  }

  #[test]
  fn test_divide_empty_by_any() {
    let empty = DistributionPlain::empty();
    let point = DistributionPlain::point(1.0, 2.0);
    let actual = distribution_division(&empty, &point, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_by_empty_fails() {
    let point = DistributionPlain::point(1.0, 2.0);
    let empty = DistributionPlain::empty();
    assert_error!(
      distribution_division(&point, &empty, MaxGridPoints::default()),
      "Cannot divide by empty distribution"
    );
  }

  #[test]
  fn test_divide_point_by_function() {
    let point = DistributionPlain::point(2.0, 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = DistributionPlain::function(t, y).unwrap();

    let actual = distribution_division(&point, &func, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::point(2.0, 2.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = DistributionPlain::function(t.clone(), y1).unwrap();
    let divisor = DistributionPlain::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();

    let expected_y = array![5.0, 5.0, 6.0, 5.0, 5.0];
    let expected = DistributionPlain::function(t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_by_zero_handled() {
    let t = array![0.0, 1.0, 2.0];
    let y1 = array![10.0, 20.0, 30.0];
    let y2 = array![2.0, 0.0, 5.0];

    let dividend = DistributionPlain::function(t.clone(), y1).unwrap();
    let divisor = DistributionPlain::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();

    let expected_y = array![5.0, 20.0 / TINY_NUMBER, 6.0];
    let expected = DistributionPlain::function(t, expected_y).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_full_overlap() {
    let range = DistributionPlain::range((1.0, 3.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y = array![1.0, 2.0, 5.0, 4.0, 3.0];
    let func = DistributionPlain::function(t, y).unwrap();

    let actual = distribution_division(&range, &func, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::Function(
      DistributionFunction::from_arrays(&array![1.0, 2.0, 3.0], array![5.0, 2.0, 2.5])
        .unwrap()
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
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
    let range = DistributionPlain::range(range_bounds, 12.0);
    let function = DistributionPlain::function(
      array![0.0, 1.0, 2.0, 3.0, 4.0],
      array![1.0, 2.0, 3.0, 4.0, 5.0],
    )
    .unwrap();

    let actual = distribution_division(&range, &function, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };
    pretty_assert_ulps_eq!(Array1::from_vec(expected_t), actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(Array1::from_vec(expected_y), actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_divide_range_by_function_endpoint_contact_returns_point() {
    let range = DistributionPlain::range((2.0, 3.0), 12.0);
    let function = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_division(&range, &function, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::point(2.0, 4.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_range_by_function_no_overlap() {
    let range = DistributionPlain::range((5.0, 6.0), 10.0);
    let t = array![0.0, 1.0, 2.0, 3.0];
    let y = array![1.0, 2.0, 4.0, 8.0];
    let func = DistributionPlain::function(t, y).unwrap();

    let actual = distribution_division(&range, &func, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::empty();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_same_grid() {
    let t = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let y2 = array![2.0, 4.0, 5.0, 8.0, 10.0];

    let dividend = DistributionPlain::function(t.clone(), y1).unwrap();
    let divisor = DistributionPlain::function(t.clone(), y2).unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();

    let expected = DistributionPlain::function(t, array![5.0, 5.0, 6.0, 5.0, 5.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_different_grids() {
    let t1 = array![0.0, 1.0, 2.0, 3.0, 4.0];
    let y1 = array![10.0, 20.0, 30.0, 40.0, 50.0];
    let dividend = DistributionPlain::function(t1.clone(), y1).unwrap();

    let t2 = array![0.0, 2.0, 4.0];
    let y2 = array![2.0, 5.0, 10.0];
    let divisor = DistributionPlain::function(t2, y2).unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();

    let expected = DistributionPlain::function(t1, array![5.0, 20.0 / 3.5, 6.0, 40.0 / 7.5, 5.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_neglog_bounds_quotient_to_hard_divisor() {
    let dividend =
      DistributionNegLog::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![0.4, 0.2, 0.0, 0.2, 0.4]).unwrap();
    let divisor = DistributionNegLog::function(array![1.0, 2.0, 3.0], array![0.1, 0.0, 0.1])
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();

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
    let dividend =
      DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![2.0, 4.0, 6.0, 8.0, 10.0]).unwrap();
    let divisor = DistributionPlain::function(array![1.5, 2.5, 3.5], array![2.0, 2.0, 2.0]).unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::function(array![1.5, 2.5, 3.5], array![2.5, 3.5, 4.5]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_soft_divisor_tail_extends_quotient_within_dividend() {
    let dividend =
      DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![10.0, 20.0, 30.0, 40.0, 50.0]).unwrap();
    let divisor = DistributionPlain::function(array![1.0, 2.0, 3.0], array![2.0, 2.0, 2.0])
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.5 }))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: 0.5 }))
      .unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();
    let expected =
      DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![4.0, 10.0, 15.0, 20.0, 20.0]).unwrap();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_endpoint_contact_returns_point() {
    let dividend = DistributionPlain::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let divisor = DistributionPlain::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_division(&dividend, &divisor, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::point(1.0, 0.6);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_divide_function_by_function_soft_left_refits_hard_right_binds() {
    let dividend =
      DistributionFunction::from_start_dx_values(2.0, 1.0, array![9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        .unwrap()
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.5 }))
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap();
    let divisor = DistributionFunction::from_start_dx_values(0.0, 1.0, Array1::from_elem(13, 2.0))
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.25 }))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_division(
      &DistributionPlain::Function(dividend),
      &DistributionPlain::Function(divisor),
      MaxGridPoints::default(),
    )
    .unwrap();

    let DistributionPlain::Function(f) = actual else {
      panic!("expected a Function quotient");
    };
    pretty_assert_ulps_eq!(0.0, f.x_min(), max_ulps = 4);
    pretty_assert_ulps_eq!(10.0, f.x_max(), max_ulps = 4);
    let BoundaryBehavior::Linear(left) = f.left_extrap() else {
      panic!("expected a refit Linear left tail, got {:?}", f.left_extrap());
    };
    assert!(
      left.slope <= 0.0,
      "refit left tail must decay, got slope {}",
      left.slope
    );
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right_end_one_ulp_outside((2.0,              3.0_f64.next_up()), (BoundaryBehavior::Hard,  BoundaryBehavior::Hard))]
  #[case::left_end_one_ulp_outside( (2.0_f64.next_down(), 3.0),            (BoundaryBehavior::Hard,  BoundaryBehavior::Hard))]
  #[case::right_end_clearly_outside((2.0,              3.5),               (BoundaryBehavior::Hard,  BoundaryBehavior::Error))]
  #[case::left_end_clearly_outside( (1.5,              3.0),               (BoundaryBehavior::Error, BoundaryBehavior::Hard))]
  #[trace]
  fn test_divide_function_by_function_coincident_hard_ends_bind_quotient(
    #[case] dividend_bounds: (f64, f64),
    #[case] expected_tails: (BoundaryBehavior, BoundaryBehavior),
  ) {
    let dividend = DistributionFunction::from_range_values(dividend_bounds, Array1::from_elem(11, 6.0))
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();
    let divisor = DistributionFunction::from_range_values((2.0, 3.0), Array1::from_elem(11, 2.0))
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_division(
      &DistributionPlain::Function(dividend),
      &DistributionPlain::Function(divisor),
      MaxGridPoints::default(),
    )
    .unwrap();

    let DistributionPlain::Function(q) = actual else {
      panic!("expected a Function quotient, got {actual:?}");
    };
    assert_eq!((2.0, 3.0), (q.x_min(), q.x_max()));
    assert_eq!(expected_tails, (q.left_extrap(), q.right_extrap()));
  }

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

    let actual = distribution_division(
      &DistributionPlain::Function(dividend),
      &DistributionPlain::Function(divisor),
      MaxGridPoints::default(),
    )
    .unwrap();

    let DistributionPlain::Function(f) = actual else {
      panic!("expected a Function quotient");
    };
    assert_eq!(BoundaryBehavior::Error, f.left_extrap());
    assert_eq!(BoundaryBehavior::Error, f.right_extrap());
  }

  #[test]
  fn test_divide_function_by_function_roundtrip_recovers_dividend_factor() {
    let grid = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    let a =
      DistributionFunction::<f64, NegLog>::from_arrays(&grid, array![4.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
        .unwrap()
        .fit_soft_tail(Side::Left, DEFAULT_TAIL_FIT_POINTS)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap();
    let b =
      DistributionFunction::<f64, NegLog>::from_arrays(&grid, array![2.0, 1.5, 1.0, 0.5, 0.0, 0.5, 1.0, 1.5, 2.0])
        .unwrap()
        .fit_soft_tail(Side::Left, DEFAULT_TAIL_FIT_POINTS)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap();

    let product = distribution_multiplication(
      &DistributionNegLog::Function(a.clone()),
      &DistributionNegLog::Function(b.clone()),
      MaxGridPoints::default(),
    )
    .unwrap();
    let quotient = distribution_division(&product, &DistributionNegLog::Function(b), MaxGridPoints::default()).unwrap();

    let DistributionNegLog::Function(q) = quotient else {
      panic!("expected a Function quotient");
    };
    pretty_assert_ulps_eq!(a.t(), q.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(a.y().clone(), q.y().clone(), max_ulps = 8);
    assert_eq!(a.left_extrap(), q.left_extrap());
    assert_eq!(a.right_extrap(), q.right_extrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::right_zero_density(Side::Right, (2.0, 0.5), f64::INFINITY,                Some(1.5))]
  #[case::right_finite(      Side::Right, (1.5, 1.5), 0.7,                          None)]
  #[case::right_divergent(   Side::Right, (0.5, 2.0), 0.7 + 1.5 * 0.1_f64.ln(),     None)]
  #[case::left_zero_density( Side::Left,  (2.0, 0.5), f64::INFINITY,                Some(1.5))]
  #[case::left_finite(       Side::Left,  (1.5, 1.5), 0.7,                          None)]
  #[case::left_divergent(    Side::Left,  (0.5, 2.0), 0.7 + 1.5 * 0.1_f64.ln(),     None)]
  #[trace]
  fn test_divide_function_by_function_shared_hard_end_takes_limit_of_approach_laws(
    #[case] side: Side,
    #[case] (dividend_exponent, divisor_exponent): (f64, f64),
    #[case] expected_end: f64,
    #[case] expected_exponent: Option<f64>,
  ) {
    let dividend = power_law_to_hard_end(side, dividend_exponent, 0.7);
    let divisor = power_law_to_hard_end(side, divisor_exponent, 0.0);

    let actual = distribution_division(
      &DistributionNegLog::Function(dividend.clone()),
      &DistributionNegLog::Function(divisor),
      MaxGridPoints::default(),
    )
    .unwrap();

    let DistributionNegLog::Function(q) = actual else {
      panic!("expected a Function quotient, got {actual:?}");
    };
    let exponent = dividend_exponent - divisor_exponent;
    let (end, mut expected) = (end_index(side), power_law_values(side, exponent, 0.7));
    expected[end] = expected_end;
    pretty_assert_ulps_eq!(dividend.t(), q.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected, q.y().clone(), max_ulps = 16);
    let expected_tail = match expected_exponent {
      Some(b) => BoundaryBehavior::HardApproach(HardApproachLaw {
        t_hard: hard_end(side),
        b,
      }),
      None => end_tail(&dividend, side),
    };
    assert_eq!(expected_tail, end_tail(&q, side));
  }

  #[test]
  fn test_divide_function_by_function_shared_hard_end_on_different_grids_stays_defined() {
    let dividend = power_law_to_hard_end(Side::Right, 2.0, 0.7)
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();
    let t = Array1::linspace(0.3, 1.0, 8);
    let y = t.mapv(|t: f64| if t < 1.0 { -(1.0 - t).ln() } else { f64::INFINITY });
    let divisor = DistributionFunction::<f64, NegLog>::from_arrays(&t, y)
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_division(
      &DistributionNegLog::Function(dividend),
      &DistributionNegLog::Function(divisor),
      MaxGridPoints::default(),
    )
    .unwrap();

    let DistributionNegLog::Function(q) = actual else {
      panic!("expected a Function quotient, got {actual:?}");
    };
    let last = q.y().len() - 1;
    assert!(
      q.y().iter().take(last).all(|v| v.is_finite()),
      "interior of the quotient must stay finite, got {:?}",
      q.y()
    );
    assert!(
      q.y()[last].is_infinite() && q.y()[last] > 0.0,
      "the shared end must have zero density, got {}",
      q.y()[last]
    );
    let BoundaryBehavior::HardApproach(law) = q.right_extrap() else {
      panic!("expected a HardApproach right tail, got {:?}", q.right_extrap());
    };
    pretty_assert_ulps_eq!(1.0, law.t_hard, max_ulps = 4);
    assert!(
      law.b > 0.0,
      "the quotient must approach zero density, got exponent {}",
      law.b
    );
  }

  mod helpers {
    use crate::__tests__::aliases::DistributionPlain;
    use crate::distribution_core::formula::DistributionFormula;
    use crate::distribution_core::function::DistributionFunction;
    use crate::policy::NegLog;
    use ndarray::{Array1, array};
    use treetime_grid::{BoundaryBehavior, HardApproachLaw, Side};
    use treetime_utils::cast::ToF64Exact;

    const N_POINTS: usize = 11;

    pub(super) fn power_law_to_hard_end(side: Side, exponent: f64, offset: f64) -> DistributionFunction<f64, NegLog> {
      let mut y = power_law_values(side, exponent, offset);
      y[end_index(side)] = f64::INFINITY;
      let law = BoundaryBehavior::HardApproach(HardApproachLaw {
        t_hard: hard_end(side),
        b: exponent,
      });
      let f = DistributionFunction::from_start_dx_values(0.0, 0.1, y).unwrap();
      match side {
        Side::Left => f
          .with_left_extrap(law)
          .unwrap()
          .with_right_extrap(BoundaryBehavior::Hard)
          .unwrap(),
        Side::Right => f
          .with_left_extrap(BoundaryBehavior::Hard)
          .unwrap()
          .with_right_extrap(law)
          .unwrap(),
      }
    }

    pub(super) fn power_law_values(side: Side, exponent: f64, offset: f64) -> Array1<f64> {
      Array1::from_shape_fn(N_POINTS, |i| {
        let dt = (0.1 * i.to_f64_exact() - hard_end(side)).abs();
        offset - exponent * dt.ln()
      })
    }

    pub(super) fn hard_end(side: Side) -> f64 {
      match side {
        Side::Left => 0.0,
        Side::Right => 1.0,
      }
    }

    pub(super) fn end_index(side: Side) -> usize {
      match side {
        Side::Left => 0,
        Side::Right => N_POINTS - 1,
      }
    }

    pub(super) fn end_tail(f: &DistributionFunction<f64, NegLog>, side: Side) -> BoundaryBehavior {
      match side {
        Side::Left => f.left_extrap(),
        Side::Right => f.right_extrap(),
      }
    }

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
