#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionPlain;
  use crate::distribution_ops::multiply::{
    HardDomain, distribution_multiplication, guarded_empty_result, hard_domains_disjoint,
  };
  use eyre::Report;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::MaxGridPoints;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, SoftTailLaw};
  use treetime_grid::{NegLog, Plain};
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};

  const SOFT: BoundaryBehavior = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.001 });

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

    let actual = distribution_multiplication(&range, &function, MaxGridPoints::default()).unwrap();
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

    let actual = distribution_multiplication(&range, &function, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::Empty;
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_range_function_endpoint_contact_returns_point() {
    let range = DistributionPlain::range((2.0, 3.0), 2.0);
    let function = DistributionPlain::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_multiplication(&range, &function, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::point(2.0, 6.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_function_function_endpoint_contact_returns_point() {
    let left = DistributionPlain::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let right = DistributionPlain::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_multiplication(&left, &right, MaxGridPoints::default()).unwrap();
    let expected = DistributionPlain::point(1.0, 15.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_function_function_uses_finer_spacing_over_intersection() {
    let coarse = DistributionPlain::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let fine = DistributionPlain::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();

    let actual = distribution_multiplication(&coarse, &fine, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    let expected_t = array![1.2, 1.7, 2.2, 2.7];
    let expected_y = array![4.4, 8.1, 12.8, 18.5];
    pretty_assert_ulps_eq!(expected_t, actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_multiply_function_function_non_overlapping_returns_empty() {
    let a = helpers::make_gaussian(0.0, 1.0, 101);
    let b = helpers::make_gaussian(20.0, 1.0, 101);

    let actual = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    assert_eq!(DistributionPlain::Empty, actual);
  }

  #[test]
  fn test_multiply_function_function_overlapping_uses_intersection() {
    let a = helpers::make_gaussian(0.0, 1.0, 101);
    let b = helpers::make_gaussian(2.0, 1.0, 101);

    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();

    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function variant");
    };

    pretty_assert_ulps_eq!(
      -3.0,
      f.x_min(),
      max_ulps = 4,
      "start of the hard support [2 - 5, 2 + 5] of N(2, 1), the inner start"
    );
    pretty_assert_ulps_eq!(
      5.0,
      f.x_max(),
      max_ulps = 4,
      "end of the hard support [0 - 5, 0 + 5] of N(0, 1), the inner end"
    );
    pretty_assert_ulps_eq!(
      1.0,
      result.likely_time().unwrap().unwrap(),
      max_ulps = 4,
      "mode of N(0, 1) N(2, 1) = N(1, 1/2) up to a constant; t = 1 = -3 + 40 * 0.1 is a node of the product grid"
    );
  }

  #[test]
  fn test_multiply_tail_c1_overlapping_no_tails() {
    let a = DistributionPlain::Function(helpers::make_function(0.0, 10.0, 101, 5.0, 2.0));
    let b = DistributionPlain::Function(helpers::make_function(3.0, 13.0, 101, 8.0, 2.0));
    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    pretty_assert_ulps_eq!(f.x_min(), 3.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
    pretty_assert_ulps_eq!(
      6.5,
      result.likely_time().unwrap().unwrap(),
      max_ulps = 4,
      "mode of N(5, 2^2) N(8, 2^2) = N(6.5, 2) up to a constant; t = 6.5 = 3 + 35 * 0.1 is a node of the product grid"
    );
  }

  #[test]
  fn test_multiply_tail_c2_overlapping_with_tails() {
    let a = DistributionPlain::Function(
      helpers::make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      helpers::make_function(3.0, 13.0, 101, 8.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    pretty_assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
  }

  #[test]
  fn test_multiply_tail_c3_disjoint_one_constant_left() {
    let a = DistributionPlain::Function(
      helpers::make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = DistributionPlain::Function(helpers::make_function(0.0, 8.0, 101, 4.0, 2.0));
    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function, got {result:?}")
    };
    pretty_assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 8.0, max_ulps = 4);
    pretty_assert_ulps_eq!(
      4.0,
      result.likely_time().unwrap().unwrap(),
      max_ulps = 4,
      "argmax on the grid of step 0.08 of exp(-(t - 4)^2 / 8 + 0.001 t), N(4, 2^2) times the exponential soft tail of a; ln p at 4.00 exceeds 4.08 by 7.2e-4 (bc -l)"
    );
  }

  #[test]
  fn test_multiply_tail_c4_disjoint_both_constant_left() {
    let leaf_msg = DistributionPlain::Function(
      helpers::make_function(2001.0, 2007.0, 61, 2004.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let subtree_msg = DistributionPlain::Function(
      helpers::make_function(1970.0, 2000.0, 301, 1990.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let result = distribution_multiplication(&leaf_msg, &subtree_msg, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function (non-empty product), got {result:?}")
    };
    pretty_assert_ulps_eq!(f.x_min(), 1970.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 2000.0, max_ulps = 4);
    pretty_assert_ulps_eq!(
      1990.0,
      result.likely_time().unwrap().unwrap(),
      max_ulps = 4,
      "argmax on the grid of step 0.1 of exp(-(t - 1990)^2 / 50 + 0.001 t), the subtree message times the exponential soft tail of the leaf message; ln p at 1990.0 exceeds 1990.1 by 1.0e-4 (bc -l)"
    );
  }

  #[test]
  fn test_multiply_tail_c5_disjoint_hard_tail() {
    let a = DistributionPlain::Function(
      helpers::make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let b = DistributionPlain::Function(helpers::make_function(0.0, 8.0, 101, 4.0, 2.0));
    let actual = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    assert_eq!(DistributionPlain::Empty, actual);
  }

  #[test]
  fn test_multiply_tail_c6_endpoint_contact_with_constant() {
    let a = DistributionPlain::Function(
      helpers::make_function(5.0, 10.0, 51, 7.5, 2.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = DistributionPlain::Function(helpers::make_function(0.0, 5.0, 51, 2.5, 2.0));
    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function (interval, not point contact), got {result:?}")
    };
    pretty_assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 5.0, max_ulps = 4);
  }

  #[test]
  fn test_multiply_tail_c7_contained_with_tails() {
    let outer = DistributionPlain::Function(
      helpers::make_function(0.0, 20.0, 201, 10.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let inner = DistributionPlain::Function(helpers::make_function(5.0, 15.0, 101, 10.0, 3.0));
    let result = distribution_multiplication(&outer, &inner, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    pretty_assert_ulps_eq!(f.x_min(), 5.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 15.0, max_ulps = 4);
  }

  #[test]
  fn test_multiply_tail_c8_mixed_constant_zero() {
    let a = DistributionPlain::Function(
      helpers::make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      helpers::make_function(3.0, 8.0, 51, 5.5, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(f) = &result else {
      panic!("Expected Function")
    };
    pretty_assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(f.x_max(), 8.0, max_ulps = 4);
    pretty_assert_ulps_eq!(
      5.3,
      result.likely_time().unwrap().unwrap(),
      max_ulps = 4,
      "argmax on the grid of step 0.1 of N(5, 2^2) N(5.5, 1.5^2) = N(5.32, 1.44) up to a constant; ln p at 5.3 exceeds 5.4 by 2.1e-3"
    );
  }

  #[test]
  fn test_multiply_tail_commutative() {
    let a = DistributionPlain::Function(
      helpers::make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = DistributionPlain::Function(helpers::make_function(0.0, 8.0, 101, 4.0, 2.0));
    let ab = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let ba = distribution_multiplication(&b, &a, MaxGridPoints::default()).unwrap();
    let (DistributionPlain::Function(fab), DistributionPlain::Function(fba)) = (&ab, &ba) else {
      panic!("Both results should be Function")
    };
    pretty_assert_ulps_eq!(fab.x_min(), fba.x_min(), max_ulps = 4);
    pretty_assert_ulps_eq!(fab.x_max(), fba.x_max(), max_ulps = 4);
    pretty_assert_ulps_eq!(fab.y(), fba.y(), max_ulps = 10);
  }

  #[test]
  fn test_multiply_tail_chained_survives_disjoint() -> Result<(), Report> {
    let [msg1, msg2, msg3, msg4] = helpers::chained_messages::<Plain>()?;
    let leaves = helpers::multiply_all(&msg1, &[&msg2, &msg3])?;
    let DistributionPlain::Function(f) = &leaves else {
      panic!("Expected Function, got {leaves:?}")
    };
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw {
        slope: -0.001 - 0.001 - 0.001
      }),
      f.left_extrap(),
      "costs add under multiplication, so the slopes of the three soft tails add"
    );

    let accum = helpers::multiply_all(&leaves, &[&msg4])?;
    let DistributionPlain::Function(f) = &accum else {
      panic!("Chained multiplication collapsed to {accum:?} despite soft left tails")
    };
    assert_eq!(
      BoundaryBehavior::Linear(SoftTailLaw {
        slope: -0.001 - 0.001 - 0.001 - 0.001
      }),
      f.left_extrap()
    );
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    pretty_assert_ulps_eq!(
      1990.1,
      accum.likely_time()?.unwrap(),
      max_ulps = 4,
      "argmax on the grid of step 0.1 of exp(-(t - 1990)^2 / 50 + 0.003 t), the subtree message times the summed soft tails of the three leaf messages; ln p at 1990.1 exceeds 1990.0 by 1.0e-4 (bc -l)"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::far_in_tail((1990.0,  7.159))]
  #[case::near_edge(  (1999.95, 7.129_15))]
  #[trace]
  fn test_multiply_tail_chained_leaf_product_is_the_product_of_the_tails(
    #[case] (t, expected): (f64, f64),
  ) -> Result<(), Report> {
    let [plain1, plain2, plain3, _] = helpers::chained_messages::<Plain>()?;
    let [neg_log1, neg_log2, neg_log3, _] = helpers::chained_messages::<NegLog>()?;
    let plain = helpers::multiply_all(&plain1, &[&plain2, &plain3])?;
    let neg_log = helpers::multiply_all(&neg_log1, &[&neg_log2, &neg_log3])?;
    pretty_assert_abs_diff_eq!(
      expected,
      neg_log.eval(t)?,
      epsilon = 1e-13,
      "cost of the product of the three tails at t < 2000: (2005 - 2000)^2 / 8 + (2004 - 2001)^2 / 4.5 + 0.001 (2001 - t) + (2006 - 2003)^2 / 4.5 + 0.001 (2003 - t) + 0.001 (2000 - t) nats. Bound: at most 30 roundings, each within one ulp (1.8e-15) of a cost below 8"
    );
    pretty_assert_abs_diff_eq!(
      expected,
      -plain.eval(t)?.ln(),
      epsilon = 1e-13,
      "the plain product is exp(-cost) of the same closed form. Bound: at most 30 roundings, each within one ulp (1.8e-15) of a cost below 8"
    );
    Ok(())
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
    let a = DistributionPlain::Function(helpers::make_function(0.0, 10.0, 101, 5.0, 2.0).with_left_extrap(a_left).unwrap());
    let b = DistributionPlain::Function(helpers::make_function(2.0, 12.0, 101, 7.0, 2.0).with_left_extrap(b_left).unwrap());

    let ab = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(fab) = &ab else {
      panic!("Expected Function, got {ab:?}")
    };
    assert_eq!(expected_left, fab.left_extrap());
    assert_eq!(BoundaryBehavior::Error, fab.right_extrap());

    let ba = distribution_multiplication(&b, &a, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(fba) = &ba else {
      panic!("Expected Function, got {ba:?}")
    };
    assert_eq!(expected_left, fba.left_extrap());
  }

  #[test]
  fn test_multiply_function_function_composes_soft_tails() {
    let a = DistributionPlain::Function(
      helpers::make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.002 }))
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      helpers::make_function(2.0, 12.0, 101, 7.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.003 }))
        .unwrap(),
    );

    let ab = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
    let DistributionPlain::Function(fab) = &ab else {
      panic!("Expected Function, got {ab:?}")
    };
    let BoundaryBehavior::Linear(law) = fab.left_extrap() else {
      panic!("Expected a composed Linear left tail, got {:?}", fab.left_extrap())
    };
    pretty_assert_ulps_eq!(-0.005, law.slope, max_ulps = 4);
    assert_eq!(BoundaryBehavior::Error, fab.right_extrap());
  }

  #[test]
  fn test_multiply_function_function_two_hard_approach_is_unsupported() {
    let law_a = HardApproachLaw { t_hard: 0.0, b: 1.0 };
    let law_b = HardApproachLaw { t_hard: 0.0, b: 2.0 };

    let a = DistributionPlain::Function(
      helpers::make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law_a))
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      helpers::make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law_b))
        .unwrap(),
    );

    assert_error!(
      distribution_multiplication(&a, &b, MaxGridPoints::default()),
      "Cannot multiply two HardApproach tails: their product is not representable by a single-parameter hard-approach law, and this composition is unreachable in the inference pipeline. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_multiply_function_function_nullary_hard_absorbs_approach_law() {
    let law = HardApproachLaw { t_hard: 0.0, b: 1.5 };

    let a = DistributionPlain::Function(
      helpers::make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law))
        .unwrap(),
    );
    let b = DistributionPlain::Function(
      helpers::make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );

    let result = distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap();
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
    let error = distribution_multiplication(&point, &func, MaxGridPoints::default())
      .unwrap_err()
      .to_string();
    assert!(error.contains("hard domains overlap"), "unexpected error: {error}");
  }

  #[test]
  fn test_multiply_disjoint_points_return_empty() {
    let a = DistributionPlain::point(1.0, 2.0);
    let b = DistributionPlain::point(5.0, 3.0);
    assert_eq!(
      DistributionPlain::Empty,
      distribution_multiplication(&a, &b, MaxGridPoints::default()).unwrap()
    );
  }

  #[test]
  fn test_guarded_empty_result_overlap_raises_internal_error() {
    let a = Some(((0.0, 10.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    let b = Some(((5.0, 15.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    assert_error!(
      guarded_empty_result::<Plain>("convolution", a, b),
      "convolution produced an empty result, but the operands' hard domains overlap: operand A Some(((0.0, 10.0), (Error, Error))), operand B Some(((5.0, 15.0), (Error, Error))). An empty result must arise only from genuinely disjoint hard domains, never from numerical collapse. This is an internal error. Please report it to developers."
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::disjoint_hard_domains( Some(((0.0, 8.0), (BoundaryBehavior::Error, BoundaryBehavior::Error))),  Some(((10.0, 20.0), (BoundaryBehavior::Error, BoundaryBehavior::Error))))]
  #[case::first_operand_empty(   None,                                                                    Some(((0.0, 10.0), (BoundaryBehavior::Error, BoundaryBehavior::Error))))]
  #[case::second_operand_empty(  Some(((0.0, 10.0), (BoundaryBehavior::Error, BoundaryBehavior::Error))), None)]
  #[trace]
  fn test_guarded_empty_result_returns_empty_without_overlap(#[case] a: Option<HardDomain>, #[case] b: Option<HardDomain>) {
    assert_eq!(
      DistributionPlain::Empty,
      guarded_empty_result::<Plain>("multiplication", a, b).unwrap()
    );
  }

  mod helpers {
    use super::SOFT;
    use crate::__tests__::aliases::DistributionPlain;
    use crate::distribution_ops::multiply::distribution_multiplication;
    use crate::{Distribution, DistributionFunction};
    use eyre::Report;
    use ndarray::Array1;
    use treetime_grid::{BoundaryBehavior, MaxGridPoints, Plain, YAxisPolicy};
    use treetime_utils::cast::ToF64Exact;

    pub(super) fn chained_messages<Y: YAxisPolicy>() -> Result<[Distribution<Y>; 4], Report> {
      Ok([
        soft_message((2000.0, 2010.0), 101, 2005.0, 2.0)?,
        soft_message((2001.0, 2008.0), 71, 2004.0, 1.5)?,
        soft_message((2003.0, 2009.0), 61, 2006.0, 1.5)?,
        soft_message((1970.0, 1999.0), 291, 1990.0, 5.0)?,
      ])
    }

    pub(super) fn multiply_all<Y: YAxisPolicy>(
      first: &Distribution<Y>,
      rest: &[&Distribution<Y>],
    ) -> Result<Distribution<Y>, Report> {
      rest.iter().try_fold(first.clone(), |accum, msg| {
        distribution_multiplication(&accum, msg, MaxGridPoints::default())
      })
    }

    fn soft_message<Y: YAxisPolicy>(
      range: (f64, f64),
      n: usize,
      peak_at: f64,
      sigma: f64,
    ) -> Result<Distribution<Y>, Report> {
      let y = Array1::linspace(range.0, range.1, n).mapv(|x| Y::from_neg_log(0.5 * ((x - peak_at) / sigma).powi(2)));
      let f = DistributionFunction::<f64, Y>::from_range_values(range, y)?
        .with_left_extrap(SOFT)?
        .with_right_extrap(BoundaryBehavior::Hard)?;
      Ok(Distribution::Function(f))
    }

    pub(super) fn make_gaussian(mu: f64, sigma: f64, n_points: usize) -> DistributionPlain {
      let x_min = mu - 5.0 * sigma;
      let x_max = mu + 5.0 * sigma;
      let dx = (x_max - x_min) / (n_points - 1).to_f64_exact();
      let y = Array1::from_shape_fn(n_points, |i| {
        let x = x_min + dx * i.to_f64_exact();
        (-0.5 * ((x - mu) / sigma).powi(2)).exp()
      });
      let f = DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y).unwrap();
      DistributionPlain::Function(f)
    }

    pub(super) fn make_function(
      x_min: f64,
      x_max: f64,
      n: usize,
      peak_at: f64,
      sigma: f64,
    ) -> DistributionFunction<f64, Plain> {
      let dx = (x_max - x_min) / (n - 1).to_f64_exact();
      let y = Array1::from_shape_fn(n, |i| {
        let x = x_min + dx * i.to_f64_exact();
        (-0.5 * ((x - peak_at) / sigma).powi(2)).exp()
      });
      DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y).unwrap()
    }
  }
}
