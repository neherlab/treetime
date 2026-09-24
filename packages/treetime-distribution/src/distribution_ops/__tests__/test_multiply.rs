#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionPlain;
  use crate::DistributionFunction;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_ops::multiply::{distribution_multiplication, guarded_empty_result, hard_domains_disjoint};
  use crate::policy::Plain;
  use approx::assert_ulps_eq;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  const SOFT: BoundaryBehavior = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.001 });

  #[test]
  fn test_multiply_formula_function_returns_function() {
    let formula = DistributionFormula::new(|t| Ok(2.0 * t), 0.0, 10.0);
    let formula_dist = DistributionPlain::Formula(formula);

    let t = array![1.0, 3.0, 5.0, 7.0, 9.0];
    let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let function_dist = DistributionPlain::function(t, y).unwrap();

    let result = distribution_multiplication(&formula_dist, &function_dist).unwrap();

    let DistributionPlain::Function(result_fn) = result else {
      panic!("Expected Function variant, got {result:?}");
    };

    let expected = array![2.0, 12.0, 30.0, 56.0, 90.0];
    assert_ulps_eq!(expected, result_fn.y(), max_ulps = 10);
  }

  #[test]
  fn test_multiply_function_formula_commutative() {
    let formula = DistributionFormula::new(|t| Ok(t * t), 0.0, 5.0);
    let formula_dist = DistributionPlain::Formula(formula);

    let t = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let function_dist = DistributionPlain::function(t, y).unwrap();

    let result_ff = distribution_multiplication(&formula_dist, &function_dist).unwrap();
    let result_fxf = distribution_multiplication(&function_dist, &formula_dist).unwrap();

    let (DistributionPlain::Function(ff_fn), DistributionPlain::Function(fxf_fn)) = (&result_ff, &result_fxf) else {
      panic!("Both results should be Function variants");
    };

    assert_ulps_eq!(ff_fn.y(), fxf_fn.y(), max_ulps = 10);
  }

  #[test]
  fn test_multiply_formula_function_uses_function_spacing_over_intersection() {
    let formula = DistributionPlain::Formula(DistributionFormula::new(|t| Ok(2.0 * t), 1.2, 2.4));
    let function =
      DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();

    let actual = distribution_multiplication(&formula, &function).unwrap();
    let DistributionPlain::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    let expected_t = array![1.2, 1.8, 2.4];
    let expected_y = array![5.28, 10.08, 16.32];
    pretty_assert_ulps_eq!(expected_t, actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::contained(        (1.0, 3.0), vec![1.0, 2.0, 3.0],                         vec![4.0, 6.0, 8.0])]
  #[case::left_partial(    (-1.0, 2.5), vec![0.0, 5.0 / 6.0, 5.0 / 3.0, 2.5],       vec![2.0, 11.0 / 3.0, 16.0 / 3.0, 7.0])]
  #[case::right_partial(    (1.5, 5.0), vec![1.5, 7.0 / 3.0, 19.0 / 6.0, 4.0],      vec![5.0, 20.0 / 3.0, 25.0 / 3.0, 10.0])]
  #[case::no_interior_knot( (1.2, 1.8), vec![1.2, 1.8],                              vec![4.4, 5.6])]
  #[trace]
  fn test_multiply_range_function_preserves_analytical_overlap(
    #[case] range_bounds: (f64, f64),
    #[case] expected_t: Vec<f64>,
    #[case] expected_y: Vec<f64>,
  ) {
    let range = DistributionPlain::range(range_bounds, 2.0);
    let function = DistributionPlain::function(
      array![0.0, 1.0, 2.0, 3.0, 4.0],
      array![1.0, 2.0, 3.0, 4.0, 5.0],
    )
    .unwrap();

    let actual = distribution_multiplication(&range, &function).unwrap();
    let DistributionPlain::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };
    pretty_assert_ulps_eq!(Array1::from_vec(expected_t), actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(Array1::from_vec(expected_y), actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_multiply_range_function_disjoint_returns_empty() {
    let range = DistributionPlain::range((5.0, 6.0), 2.0);
    let function = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_multiplication(&range, &function).unwrap();
    let expected = DistributionPlain::Empty;
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_range_function_endpoint_contact_returns_point() {
    let range = DistributionPlain::range((2.0, 3.0), 2.0);
    let function = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_multiplication(&range, &function).unwrap();
    let expected = DistributionPlain::point(2.0, 6.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_function_function_endpoint_contact_returns_point() {
    let left = DistributionPlain::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let right = DistributionPlain::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_multiplication(&left, &right).unwrap();
    let expected = DistributionPlain::point(1.0, 15.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_function_function_uses_finer_spacing_over_intersection() {
    let coarse = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let fine = DistributionPlain::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();

    let actual = distribution_multiplication(&coarse, &fine).unwrap();
    let DistributionPlain::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    let expected_t = array![1.2, 1.7, 2.2, 2.7];
    let expected_y = array![4.4, 8.1, 12.8, 18.5];
    pretty_assert_ulps_eq!(expected_t, actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  fn make_gaussian(mu: f64, sigma: f64, n_points: usize) -> DistributionPlain {
    let x_min = mu - 5.0 * sigma;
    let x_max = mu + 5.0 * sigma;
    let dx = (x_max - x_min) / (n_points - 1) as f64;
    let y = Array1::from_shape_fn(n_points, |i| {
      let x = x_min + dx * (i as f64);
      (-0.5 * ((x - mu) / sigma).powi(2)).exp()
    });
    let f = DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y).unwrap();
    DistributionPlain::Function(f)
  }

  #[test]
  fn test_multiply_function_function_non_overlapping_returns_empty() {
    let a = make_gaussian(0.0, 1.0, 101);
    let b = make_gaussian(20.0, 1.0, 101);

    let result = distribution_multiplication(&a, &b).unwrap();
    assert!(matches!(result, DistributionPlain::Empty));
  }

  #[test]
  fn test_multiply_function_function_overlapping_uses_intersection() {
    let a = make_gaussian(0.0, 1.0, 101);
    let b = make_gaussian(2.0, 1.0, 101);

    let result = distribution_multiplication(&a, &b).unwrap();

    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function variant");
    };

    let overlap_min = (0.0_f64 - 5.0).max(2.0 - 5.0);
    let overlap_max = (0.0_f64 + 5.0).min(2.0 + 5.0);
    assert!(overlap_min < overlap_max, "Distributions must overlap");

    assert!(f.x_min() >= overlap_min - 0.1);
    assert!(f.x_max() <= overlap_max + 0.1);

    let likely = result.likely_time().unwrap().unwrap();
    assert!((likely - 1.0).abs() < 0.2, "Product peak at {likely}, expected ~1.0");
  }

  fn make_function(x_min: f64, x_max: f64, n: usize, peak_at: f64, sigma: f64) -> DistributionFunction<f64, Plain> {
    let dx = (x_max - x_min) / (n - 1) as f64;
    let y = Array1::from_shape_fn(n, |i| {
      let x = x_min + dx * (i as f64);
      (-0.5 * ((x - peak_at) / sigma).powi(2)).exp()
    });
    DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y).unwrap()
  }

  #[test]
  fn test_multiply_tail_c1_overlapping_no_tails() {
    let a = DistributionPlain::Function(make_function(0.0, 10.0, 101, 5.0, 2.0));
    let b = DistributionPlain::Function(make_function(3.0, 13.0, 101, 8.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    assert_ulps_eq!(f.x_min(), 3.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
    let peak = result.likely_time().unwrap().unwrap();
    assert!(peak > 5.5 && peak < 7.5, "Peak at {peak}, expected ~6.5");
  }

  #[test]
  fn test_multiply_tail_c2_overlapping_with_tails() {
    let a = DistributionPlain::Function(make_function(0.0, 10.0, 101, 5.0, 2.0).with_left_extrap(SOFT).unwrap());
    let b = DistributionPlain::Function(make_function(3.0, 13.0, 101, 8.0, 2.0).with_left_extrap(SOFT).unwrap());
    let result = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
  }

  #[test]
  fn test_multiply_tail_c3_disjoint_one_constant_left() {
    let a = DistributionPlain::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = DistributionPlain::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function, got {result:?}")
    };
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 8.0, max_ulps = 4);
    assert!(result.likely_time().unwrap().is_some());
  }

  #[test]
  fn test_multiply_tail_c4_disjoint_both_constant_left() {
    let leaf_msg = DistributionPlain::Function(
      make_function(2001.0, 2007.0, 61, 2004.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let subtree_msg = DistributionPlain::Function(
      make_function(1970.0, 2000.0, 301, 1990.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let result = distribution_multiplication(&leaf_msg, &subtree_msg).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function (non-empty product), got {result:?}")
    };
    assert_ulps_eq!(f.x_min(), 1970.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 2000.0, max_ulps = 4);
    let peak = result.likely_time().unwrap().unwrap();
    assert!(peak > 1985.0 && peak < 1995.0, "Peak at {peak}, expected near 1990");
  }

  #[test]
  fn test_multiply_tail_c5_disjoint_hard_tail() {
    let a = DistributionPlain::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let b = DistributionPlain::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    assert!(
      matches!(result, DistributionPlain::Empty),
      "Hard tail should not prevent Empty"
    );
  }

  #[test]
  fn test_multiply_tail_c6_endpoint_contact_with_constant() {
    let a = DistributionPlain::Function(make_function(5.0, 10.0, 51, 7.5, 2.0).with_left_extrap(SOFT).unwrap());
    let b = DistributionPlain::Function(make_function(0.0, 5.0, 51, 2.5, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function (interval, not point contact), got {result:?}")
    };
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 5.0, max_ulps = 4);
  }

  #[test]
  fn test_multiply_tail_c7_contained_with_tails() {
    let outer = DistributionPlain::Function(make_function(0.0, 20.0, 201, 10.0, 5.0).with_left_extrap(SOFT).unwrap());
    let inner = DistributionPlain::Function(make_function(5.0, 15.0, 101, 10.0, 3.0));
    let result = distribution_multiplication(&outer, &inner).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    assert_ulps_eq!(f.x_min(), 5.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 15.0, max_ulps = 4);
  }

  #[test]
  fn test_multiply_tail_c8_mixed_constant_zero() {
    let a = DistributionPlain::Function(
      make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      make_function(3.0, 8.0, 51, 5.5, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let result = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 8.0, max_ulps = 4);
    assert!(result.likely_time().unwrap().is_some());
  }

  #[test]
  fn test_multiply_tail_commutative() {
    let a = DistributionPlain::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = DistributionPlain::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let ab = distribution_multiplication(&a, &b).unwrap();
    let ba = distribution_multiplication(&b, &a).unwrap();
    let (DistributionPlain::Function(fab), DistributionPlain::Function(fba)) = (&ab, &ba) else {
      panic!("Both results should be Function")
    };
    assert_ulps_eq!(fab.x_min(), fba.x_min(), max_ulps = 4);
    assert_ulps_eq!(fab.x_max(), fba.x_max(), max_ulps = 4);
    assert_ulps_eq!(fab.y(), fba.y(), max_ulps = 10);
  }

  #[test]
  fn test_multiply_tail_chained_survives_disjoint() {
    let msg1 = DistributionPlain::Function(
      make_function(2000.0, 2010.0, 101, 2005.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg2 = DistributionPlain::Function(
      make_function(2001.0, 2008.0, 71, 2004.0, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg3 = DistributionPlain::Function(
      make_function(2003.0, 2009.0, 61, 2006.0, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg4 = DistributionPlain::Function(
      make_function(1970.0, 1999.0, 291, 1990.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );

    let mut accum = msg1;
    for msg in [&msg2, &msg3, &msg4] {
      accum = distribution_multiplication(&accum, msg).unwrap();
    }

    let DistributionPlain::Function(f) = &accum else {
      panic!("Chained multiplication collapsed to {accum:?} despite soft left tails")
    };
    assert!(matches!(f.left_extrap(), BoundaryBehavior::Linear(_)));
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    assert!(
      accum.likely_time().unwrap().is_some(),
      "Accumulated result must have a likely_time"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::soft_hard(  SOFT,                    BoundaryBehavior::Hard,  BoundaryBehavior::Hard)]
  #[case::soft_error( SOFT,                    BoundaryBehavior::Error, BoundaryBehavior::Error)]
  #[case::hard_hard(  BoundaryBehavior::Hard,  BoundaryBehavior::Hard,  BoundaryBehavior::Hard)]
  #[case::hard_error( BoundaryBehavior::Hard,  BoundaryBehavior::Error, BoundaryBehavior::Error)]
  #[case::error_error(BoundaryBehavior::Error, BoundaryBehavior::Error, BoundaryBehavior::Error)]
  #[trace]
  fn test_multiply_function_function_composes_result_tails(
    #[case] a_left: BoundaryBehavior,
    #[case] b_left: BoundaryBehavior,
    #[case] expected_left: BoundaryBehavior,
  ) {
    let a = DistributionPlain::Function(make_function(0.0, 10.0, 101, 5.0, 2.0).with_left_extrap(a_left).unwrap());
    let b = DistributionPlain::Function(make_function(2.0, 12.0, 101, 7.0, 2.0).with_left_extrap(b_left).unwrap());

    let ab = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(fab) = &ab else {
      panic!("Expected Function, got {ab:?}")
    };
    assert_eq!(expected_left, fab.left_extrap());
    assert_eq!(BoundaryBehavior::Error, fab.right_extrap());

    let ba = distribution_multiplication(&b, &a).unwrap();
    let DistributionPlain::Function(fba) = &ba else {
      panic!("Expected Function, got {ba:?}")
    };
    assert_eq!(expected_left, fba.left_extrap());
  }

  #[test]
  fn test_multiply_function_function_composes_soft_tails() {
    let a = DistributionPlain::Function(
      make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.002 }))
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      make_function(2.0, 12.0, 101, 7.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.003 }))
        .unwrap(),
    );

    let ab = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(fab) = &ab else {
      panic!("Expected Function, got {ab:?}")
    };
    let BoundaryBehavior::Linear(law) = fab.left_extrap() else {
      panic!("Expected a composed Linear left tail, got {:?}", fab.left_extrap())
    };
    assert_ulps_eq!(-0.005, law.slope, max_ulps = 4);
    assert_eq!(BoundaryBehavior::Error, fab.right_extrap());
  }

  #[test]
  fn test_multiply_function_function_two_hard_approach_is_unsupported() {
    let law_a = HardApproachLaw { t_hard: 0.0, b: 1.0 };
    let law_b = HardApproachLaw { t_hard: 0.0, b: 2.0 };

    let a = DistributionPlain::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law_a))
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law_b))
        .unwrap(),
    );

    assert_error!(
      distribution_multiplication(&a, &b),
      "Cannot multiply two HardApproach tails: their product is not representable by a single-parameter hard-approach law, and this composition is unreachable in the inference pipeline. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_multiply_function_function_nullary_hard_absorbs_approach_law() {
    let law = HardApproachLaw { t_hard: 0.0, b: 1.5 };

    let a = DistributionPlain::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law))
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );

    let result = distribution_multiplication(&a, &b).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    assert_eq!(BoundaryBehavior::Hard, f.left_extrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::hard_gap_disjoint(          (10.0, 20.0, BoundaryBehavior::Error, BoundaryBehavior::Error), ( 0.0,  8.0, BoundaryBehavior::Error, BoundaryBehavior::Error), true)]
  #[case::overlapping(                ( 0.0, 10.0, BoundaryBehavior::Error, BoundaryBehavior::Error), ( 5.0, 15.0, BoundaryBehavior::Error, BoundaryBehavior::Error), false)]
  #[case::endpoint_contact(           ( 0.0, 10.0, BoundaryBehavior::Error, BoundaryBehavior::Error), (10.0, 20.0, BoundaryBehavior::Error, BoundaryBehavior::Error), false)]
  #[case::soft_left_bridges_gap(      (10.0, 20.0, SOFT, BoundaryBehavior::Error),                   ( 0.0,  8.0, BoundaryBehavior::Error, BoundaryBehavior::Error), false)]
  #[case::facing_bounds_hard_far_soft((10.0, 20.0, BoundaryBehavior::Error, SOFT),                   ( 0.0,  8.0, SOFT, BoundaryBehavior::Error),                   true)]
  #[case::backward_messages_never(    (2001.0, 2007.0, SOFT, BoundaryBehavior::Hard),                (1970.0, 2000.0, SOFT, BoundaryBehavior::Hard),                false)]
  #[trace]
  fn test_multiply_hard_domains_disjoint(
    #[case] (a_lo, a_hi, a_left, a_right): (f64, f64, BoundaryBehavior, BoundaryBehavior),
    #[case] (b_lo, b_hi, b_left, b_right): (f64, f64, BoundaryBehavior, BoundaryBehavior),
    #[case] expected: bool,
  ) {
    let actual = hard_domains_disjoint((a_lo, a_hi), (a_left, a_right), (b_lo, b_hi), (b_left, b_right));
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_point_on_function_zero_raises_internal_error() {
    let point = DistributionPlain::point(1.0, 5.0);
    let func = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 0.0, 1.0]).unwrap();
    let error = distribution_multiplication(&point, &func).unwrap_err().to_string();
    assert!(error.contains("hard domains overlap"), "unexpected error: {error}");
  }

  #[test]
  fn test_multiply_disjoint_points_return_empty() {
    let a = DistributionPlain::point(1.0, 2.0);
    let b = DistributionPlain::point(5.0, 3.0);
    assert_eq!(DistributionPlain::Empty, distribution_multiplication(&a, &b).unwrap());
  }

  #[rstest]
  #[case::multiplication("multiplication")]
  #[case::convolution("convolution")]
  #[case::division("division")]
  fn test_guarded_empty_result_overlap_raises_internal_error(#[case] operation: &str) {
    let a = Some(((0.0, 10.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    let b = Some(((5.0, 15.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    let error = guarded_empty_result::<Plain>(operation, a, b).unwrap_err().to_string();
    assert!(error.contains(operation), "message should name the operation: {error}");
    assert!(error.contains("hard domains overlap"), "unexpected error: {error}");
  }

  #[rstest]
  #[case::multiplication("multiplication")]
  #[case::convolution("convolution")]
  #[case::division("division")]
  fn test_guarded_empty_result_disjoint_returns_empty(#[case] operation: &str) {
    let a = Some(((0.0, 8.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    let b = Some(((10.0, 20.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    assert_eq!(
      DistributionPlain::Empty,
      guarded_empty_result::<Plain>(operation, a, b).unwrap()
    );
  }

  #[test]
  fn test_guarded_empty_result_empty_operand_returns_empty() {
    let present = Some(((0.0, 10.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    assert_eq!(
      DistributionPlain::Empty,
      guarded_empty_result::<Plain>("multiplication", None, present).unwrap()
    );
    assert_eq!(
      DistributionPlain::Empty,
      guarded_empty_result::<Plain>("multiplication", present, None).unwrap()
    );
  }
}
