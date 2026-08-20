#[cfg(test)]
mod tests {
  use crate::DistributionFunction;
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_ops::multiply::{distribution_multiplication, guarded_empty_result, hard_domains_disjoint};
  use crate::policy::Plain;
  use approx::assert_abs_diff_eq;
  use approx::assert_ulps_eq;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::{BoundaryBehavior, HardApproachLaw, SoftTailLaw};
  use treetime_utils::{assert_error, pretty_assert_ulps_eq};

  /// A gentle soft `Linear` left tail, the fixture stand-in for a backward-pass parent-older tail.
  /// `Linear` is the only soft `BoundaryBehavior`, so it is what bridges disjoint grids under
  /// multiplication. The slope is tiny and negative so the extrapolated Plain ordinates stay positive
  /// across the wide gaps these fixtures span (a steeper tail would underflow to a non-positive
  /// amplitude before reaching the far operand).
  const SOFT: BoundaryBehavior = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.001 });

  /// Formula * Function returns a Function with correct pointwise products.
  ///
  /// Regression test for the fix in commit b3d97b34: `multiply_formula_function`
  /// must evaluate the Formula on the Function's grid and return a concrete
  /// Function, not a lazy Formula that panics on grid-based operations.
  #[test]
  fn test_multiply_formula_function_returns_function() {
    // Formula: f(t) = 2*t over [0, 10]
    let formula = DistributionFormula::new(|t| Ok(2.0 * t), 0.0, 10.0);
    let formula_dist = Distribution::Formula(formula);

    // Function on a 5-point grid [1, 3, 5, 7, 9] with values [1, 2, 3, 4, 5]
    let t = array![1.0, 3.0, 5.0, 7.0, 9.0];
    let y = array![1.0, 2.0, 3.0, 4.0, 5.0];
    let function_dist = Distribution::function(t, y).unwrap();

    let result = distribution_multiplication(&formula_dist, &function_dist).unwrap();

    // Result must be a Function, not a Formula
    let Distribution::Function(result_fn) = result else {
      panic!("Expected Function variant, got {result:?}");
    };

    // The grid is resampled to 5 points over [1, 9].
    // At each grid point t_i, the product is formula(t_i) * function(t_i) = 2*t_i * interp(t_i).
    // Grid: t_i = 1 + (9-1)*i/4 = 1, 3, 5, 7, 9
    // formula(1)=2, formula(3)=6, formula(5)=10, formula(7)=14, formula(9)=18
    // function values: 1, 2, 3, 4, 5
    // products: 2, 12, 30, 56, 90
    let expected = array![2.0, 12.0, 30.0, 56.0, 90.0];
    assert_ulps_eq!(expected, result_fn.y(), max_ulps = 10);
  }

  /// Commutativity: Function * Formula gives the same result as Formula * Function.
  #[test]
  fn test_multiply_function_formula_commutative() {
    let formula = DistributionFormula::new(|t| Ok(t * t), 0.0, 5.0);
    let formula_dist = Distribution::Formula(formula);

    let t = array![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let y = array![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let function_dist = Distribution::function(t, y).unwrap();

    let result_ff = distribution_multiplication(&formula_dist, &function_dist).unwrap();
    let result_fxf = distribution_multiplication(&function_dist, &formula_dist).unwrap();

    let (Distribution::Function(ff_fn), Distribution::Function(fxf_fn)) = (&result_ff, &result_fxf) else {
      panic!("Both results should be Function variants");
    };

    assert_ulps_eq!(ff_fn.y(), fxf_fn.y(), max_ulps = 10);
  }

  #[test]
  fn test_multiply_formula_function_uses_function_spacing_over_intersection() {
    let formula = Distribution::Formula(DistributionFormula::new(|t| Ok(2.0 * t), 1.2, 2.4));
    let function = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();

    let actual = distribution_multiplication(&formula, &function).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    // Oracle: Richard Neher's grid contract in commit 542ac860c7cfa4bab6764aee1d1b3810a09eb54f:
    // round(overlap_width / function.dx()) + 1 points over the exact intersection.
    let expected_t = array![1.2, 2.4];
    let expected_y = array![5.28, 16.32];
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
    let range = Distribution::range(range_bounds, 2.0);
    let function = Distribution::function(
      array![0.0, 1.0, 2.0, 3.0, 4.0],
      array![1.0, 2.0, 3.0, 4.0, 5.0],
    )
    .unwrap();

    let actual = distribution_multiplication(&range, &function).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };
    pretty_assert_ulps_eq!(Array1::from_vec(expected_t), actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(Array1::from_vec(expected_y), actual.y(), max_ulps = 4);
  }

  #[test]
  fn test_multiply_range_function_disjoint_returns_empty() {
    let range = Distribution::range((5.0, 6.0), 2.0);
    let function = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_multiplication(&range, &function).unwrap();
    let expected = Distribution::Empty;
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_range_function_endpoint_contact_returns_point() {
    let range = Distribution::range((2.0, 3.0), 2.0);
    let function = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 2.0, 3.0]).unwrap();

    let actual = distribution_multiplication(&range, &function).unwrap();
    // Oracle: v0 `Distribution.multiply()` converts one surviving endpoint knot to a delta.
    let expected = Distribution::point(2.0, 6.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_function_function_endpoint_contact_returns_point() {
    let left = Distribution::function(array![0.0, 1.0], array![2.0, 3.0]).unwrap();
    let right = Distribution::function(array![1.0, 2.0], array![5.0, 7.0]).unwrap();

    let actual = distribution_multiplication(&left, &right).unwrap();
    // Oracle: v0 `Distribution.multiply()` converts one surviving endpoint knot to a delta.
    let expected = Distribution::point(1.0, 15.0);
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_multiply_function_function_uses_finer_spacing_over_intersection() {
    let coarse = Distribution::function(array![0.0, 1.0, 2.0, 3.0, 4.0], array![1.0, 2.0, 3.0, 4.0, 5.0]).unwrap();
    let fine = Distribution::function(array![1.2, 1.7, 2.2, 2.7], array![2.0, 3.0, 4.0, 5.0]).unwrap();

    let actual = distribution_multiplication(&coarse, &fine).unwrap();
    let Distribution::Function(actual) = actual else {
      panic!("Expected Function variant, got {actual:?}");
    };

    // Oracle: Richard Neher's grid contract in commit 542ac860c7cfa4bab6764aee1d1b3810a09eb54f:
    // round(overlap_width / min(dx)) + 1 points over the exact intersection.
    let expected_t = array![1.2, 1.7, 2.2, 2.7];
    let expected_y = array![4.4, 8.1, 12.8, 18.5];
    pretty_assert_ulps_eq!(expected_t, actual.t(), max_ulps = 4);
    pretty_assert_ulps_eq!(expected_y, actual.y(), max_ulps = 4);
  }

  fn make_gaussian(mu: f64, sigma: f64, n_points: usize) -> Distribution {
    let x_min = mu - 5.0 * sigma;
    let x_max = mu + 5.0 * sigma;
    let dx = (x_max - x_min) / (n_points - 1) as f64;
    let y = Array1::from_shape_fn(n_points, |i| {
      let x = x_min + dx * (i as f64);
      (-0.5 * ((x - mu) / sigma).powi(2)).exp()
    });
    let f = DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y).unwrap();
    Distribution::Function(f)
  }

  /// Non-overlapping Function * Function with default (Error) tails returns Empty.
  #[test]
  fn test_multiply_function_function_non_overlapping_returns_empty() {
    let a = make_gaussian(0.0, 1.0, 101);
    let b = make_gaussian(20.0, 1.0, 101);

    let result = distribution_multiplication(&a, &b).unwrap();
    assert!(matches!(result, Distribution::Empty));
  }

  /// Overlapping Function × Function uses intersection range.
  #[test]
  fn test_multiply_function_function_overlapping_uses_intersection() {
    // Two Gaussians 2 sigma apart: significant overlap
    let a = make_gaussian(0.0, 1.0, 101);
    let b = make_gaussian(2.0, 1.0, 101);

    let result = distribution_multiplication(&a, &b).unwrap();

    let Distribution::Function(f) = &result else {
      panic!("Expected Function variant");
    };

    // Intersection range: [2-5, 0+5] = [-3, 5]
    let overlap_min = (0.0_f64 - 5.0).max(2.0 - 5.0);
    let overlap_max = (0.0_f64 + 5.0).min(2.0 + 5.0);
    assert!(overlap_min < overlap_max, "Distributions must overlap");

    // Result range should be the intersection, not the union
    assert!(f.x_min() >= overlap_min - 0.1);
    assert!(f.x_max() <= overlap_max + 0.1);

    // Peak should be near the analytical mean of the product: mu* = 1.0
    let likely = result.likely_time().unwrap();
    assert!((likely - 1.0).abs() < 0.2, "Product peak at {likely}, expected ~1.0");
  }

  /// Chain multiplication with normalization does not underflow.
  ///
  /// Simulates the backward pass pattern: each step multiplies the
  /// accumulated product (normalized to max=1.0) with a new message.
  /// Without normalization, 50 multiplications would underflow.
  #[test]
  fn test_multiply_chain_with_normalize_no_underflow() {
    let mut accum = make_gaussian(0.0, 1.0, 201);

    for i in 0..50 {
      let msg = make_gaussian(0.0 + i as f64 * 0.01, 1.0, 201);
      accum = distribution_multiplication(&accum, &msg).unwrap().normalize();

      assert!(
        !matches!(accum, Distribution::Empty),
        "Chain multiplication underflowed at step {i}"
      );
    }

    let likely = accum.likely_time().expect("Final distribution must have likely_time");
    // After 50 steps with shifts of 0.01, the accumulated peak drifts slightly
    assert!(likely.abs() < 1.0, "Peak at {likely}, expected near 0");
  }

  fn make_function(x_min: f64, x_max: f64, n: usize, peak_at: f64, sigma: f64) -> DistributionFunction<f64, Plain> {
    let dx = (x_max - x_min) / (n - 1) as f64;
    let y = Array1::from_shape_fn(n, |i| {
      let x = x_min + dx * (i as f64);
      (-0.5 * ((x - peak_at) / sigma).powi(2)).exp()
    });
    DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y).unwrap()
  }

  /// C1: overlapping grids, no tails -- unchanged from pre-tail behavior.
  #[test]
  fn test_multiply_tail_c1_overlapping_no_tails() {
    let a = Distribution::Function(make_function(0.0, 10.0, 101, 5.0, 2.0));
    let b = Distribution::Function(make_function(3.0, 13.0, 101, 8.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // Intersection: [3.0, 10.0]
    assert_ulps_eq!(f.x_min(), 3.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
    let peak = result.likely_time().unwrap();
    // Product of two Gaussians (mu=5, mu=8) peaks at the weighted mean ~6.5
    assert!(peak > 5.5 && peak < 7.5, "Peak at {peak}, expected ~6.5");
  }

  /// C2: overlapping grids with soft tails -- intersection unchanged (tails only extend).
  #[test]
  fn test_multiply_tail_c2_overlapping_with_tails() {
    let a = Distribution::Function(make_function(0.0, 10.0, 101, 5.0, 2.0).with_left_extrap(SOFT).unwrap());
    let b = Distribution::Function(make_function(3.0, 13.0, 101, 8.0, 2.0).with_left_extrap(SOFT).unwrap());
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // A's soft left tail extends to min(0, 3) = 0 (no change). B's extends to min(3, 0) = 0.
    // Intersection of [0, 10] and [0, 13] = [0, 10].
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
  }

  /// C3: disjoint grids, one soft left tail -- extends to produce non-empty result.
  #[test]
  fn test_multiply_tail_c3_disjoint_one_constant_left() {
    // A: grid [10, 20], soft left tail
    // B: grid [0, 8], default Error tails
    let a = Distribution::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = Distribution::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function, got {result:?}")
    };
    // A extends left to min(10, 0) = 0. B stays [0, 8].
    // Intersection of [0, 20] and [0, 8] = [0, 8].
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 8.0, max_ulps = 4);
    assert!(result.likely_time().is_some());
  }

  /// C4: disjoint grids, both soft left -- backward-pass pattern.
  #[test]
  fn test_multiply_tail_c4_disjoint_both_constant_left() {
    // Two backward-pass messages with disjoint finite grids but overlapping via their soft left
    // tails: leaf message [2001, 2007] and subtree message [1970, 2000].
    let leaf_msg = Distribution::Function(
      make_function(2001.0, 2007.0, 61, 2004.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let subtree_msg = Distribution::Function(
      make_function(1970.0, 2000.0, 301, 1990.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let result = distribution_multiplication(&leaf_msg, &subtree_msg).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function (non-empty product), got {result:?}")
    };
    // Leaf extends left to min(2001, 1970) = 1970. Subtree extends left to min(1970, 2001) = 1970.
    // Right sides: both Hard, no extension.
    // Intersection of [1970, 2007] and [1970, 2000] = [1970, 2000].
    assert_ulps_eq!(f.x_min(), 1970.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 2000.0, max_ulps = 4);
    // Product should be dominated by the subtree's shape (the leaf's soft tail varies slowly here).
    let peak = result.likely_time().unwrap();
    assert!(peak > 1985.0 && peak < 1995.0, "Peak at {peak}, expected near 1990");
  }

  /// C5: disjoint grids, Hard tail -- stays Empty (zero outside support, intersection correct).
  #[test]
  fn test_multiply_tail_c5_disjoint_hard_tail() {
    let a = Distribution::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let b = Distribution::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    assert!(
      matches!(result, Distribution::Empty),
      "Hard tail should not prevent Empty"
    );
  }

  /// C6: endpoint contact with a soft tail -- extends past contact to produce interval.
  #[test]
  fn test_multiply_tail_c6_endpoint_contact_with_constant() {
    // Without tails: endpoint contact -> Point. With a soft left tail on A: extends past.
    let a = Distribution::Function(make_function(5.0, 10.0, 51, 7.5, 2.0).with_left_extrap(SOFT).unwrap());
    let b = Distribution::Function(make_function(0.0, 5.0, 51, 2.5, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function (interval, not point contact), got {result:?}")
    };
    // A extends left to min(5, 0) = 0. B stays [0, 5].
    // Intersection of [0, 10] and [0, 5] = [0, 5].
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 5.0, max_ulps = 4);
  }

  /// C7: one contained in other with tails -- intersection unchanged.
  #[test]
  fn test_multiply_tail_c7_contained_with_tails() {
    let outer = Distribution::Function(make_function(0.0, 20.0, 201, 10.0, 5.0).with_left_extrap(SOFT).unwrap());
    let inner = Distribution::Function(make_function(5.0, 15.0, 101, 10.0, 3.0));
    let result = distribution_multiplication(&outer, &inner).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // Inner is fully contained -- tails don't change the intersection.
    assert_ulps_eq!(f.x_min(), 5.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 15.0, max_ulps = 4);
  }

  /// C8: mixed soft left + Hard right -- backward message pattern.
  #[test]
  fn test_multiply_tail_c8_mixed_constant_zero() {
    // Both messages have soft left (parent could be older) and Hard right (hard upper bound).
    // Overlapping case -- verify tails don't corrupt the result.
    let a = Distribution::Function(
      make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let b = Distribution::Function(
      make_function(3.0, 8.0, 51, 5.5, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // Both extend left to min of both = 0. Right: both Hard, no extension.
    // Intersection of [0, 10] and [0, 8] = [0, 8].
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 8.0, max_ulps = 4);
    assert!(result.likely_time().is_some());
  }

  /// Commutativity holds with tail-extended multiplication.
  #[test]
  fn test_multiply_tail_commutative() {
    let a = Distribution::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(SOFT)
        .unwrap(),
    );
    let b = Distribution::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let ab = distribution_multiplication(&a, &b).unwrap();
    let ba = distribution_multiplication(&b, &a).unwrap();
    let (Distribution::Function(fab), Distribution::Function(fba)) = (&ab, &ba) else {
      panic!("Both results should be Function")
    };
    assert_ulps_eq!(fab.x_min(), fba.x_min(), max_ulps = 4);
    assert_ulps_eq!(fab.x_max(), fba.x_max(), max_ulps = 4);
    assert_ulps_eq!(fab.y(), fba.y(), max_ulps = 10);
  }

  /// Chained multiply + normalize (no manual re-apply) survives a disjoint late child.
  ///
  /// Simulates the backward pass: three overlapping messages produce a narrow accumulated
  /// result, then a fourth message with a disjoint finite grid but overlapping via its soft left tail. Multiplication
  /// left tail. Multiplication composes the result tails from its operands and normalize()
  /// preserves them, so the accumulator keeps its soft left tail across steps without any
  /// manual re-application. The fourth multiplication therefore extends leftward and stays
  /// non-empty.
  #[test]
  fn test_multiply_tail_chained_normalize_preserves_tails_survives_disjoint() {
    let msg1 = Distribution::Function(
      make_function(2000.0, 2010.0, 101, 2005.0, 2.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg2 = Distribution::Function(
      make_function(2001.0, 2008.0, 71, 2004.0, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg3 = Distribution::Function(
      make_function(2003.0, 2009.0, 61, 2006.0, 1.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    // Fourth message: disjoint from the narrowed accumulator but reachable via soft left tail.
    let msg4 = Distribution::Function(
      make_function(1970.0, 1999.0, 291, 1990.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );

    // Backward-pass accumulation: multiply + normalize only, exactly as the inference pass runs.
    let mut accum = msg1;
    for msg in [&msg2, &msg3, &msg4] {
      accum = distribution_multiplication(&accum, msg).unwrap().normalize();
    }

    let Distribution::Function(f) = &accum else {
      panic!("Chained multiplication collapsed to {accum:?} despite soft left tails")
    };
    // The composed tails survive normalization across every step.
    assert!(matches!(f.left_extrap(), BoundaryBehavior::Linear(_)));
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    assert!(
      accum.likely_time().is_some(),
      "Accumulated result must have a likely_time"
    );
  }

  /// The specific narrow-then-disjoint scenario that used to collapse now survives.
  ///
  /// Two recent messages narrow the accumulator to ~(2024.5, 2025.5). A much older message at
  /// (1970, 1999) has a disjoint finite grid. Because normalize() preserves the accumulator's
  /// soft left tail, the final multiplication extends leftward to reach the old message and
  /// produces a non-empty result. Under the old behavior, where normalize reset tails to Error,
  /// this product collapsed to Empty.
  #[test]
  fn test_multiply_tail_chained_normalize_no_reapply_survives() {
    let msg_recent_1 = Distribution::Function(
      make_function(2024.0, 2026.0, 21, 2025.0, 0.5)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg_recent_2 = Distribution::Function(
      make_function(2024.5, 2025.5, 11, 2025.0, 0.3)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );
    let msg_old = Distribution::Function(
      make_function(1970.0, 1999.0, 291, 1990.0, 5.0)
        .with_left_extrap(SOFT)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );

    // First two messages overlap and narrow the accumulator; the product keeps a soft left tail.
    let product = distribution_multiplication(&msg_recent_1, &msg_recent_2).unwrap();
    assert!(!matches!(product, Distribution::Empty));

    // Normalize preserves the soft left tail.
    let normalized = product.normalize();
    let Distribution::Function(f) = &normalized else {
      panic!("Expected Function, got {normalized:?}")
    };
    assert!(matches!(f.left_extrap(), BoundaryBehavior::Linear(_)));

    // The old message's disjoint grid is reachable because the accumulator extends leftward.
    let final_result = distribution_multiplication(&normalized, &msg_old).unwrap();
    assert!(
      !matches!(final_result, Distribution::Empty),
      "soft left tail must survive normalize() so the accumulator overlaps the old message"
    );
  }

  /// normalize() preserves the per-side tail policies of a Function distribution.
  #[test]
  fn test_multiply_normalize_preserves_tails() {
    let f = make_function(2000.0, 2010.0, 101, 2005.0, 2.0)
      .with_left_extrap(SOFT)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();
    let normalized = Distribution::Function(f).normalize();
    let Distribution::Function(f) = &normalized else {
      panic!("Expected Function, got {normalized:?}")
    };
    assert!(matches!(f.left_extrap(), BoundaryBehavior::Linear(_)));
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
    // Peak scaled to 1.0, tails unchanged.
    assert_ulps_eq!(1.0, f.y().iter().copied().fold(f64::MIN, f64::max), max_ulps = 4);
  }

  /// Function * Function composes the per-side result tail from the operands' tails, with precedence
  /// soft (`Linear`) < `Hard` < `Error` (the more restrictive tail wins). A hard bound dominates a
  /// soft one and an undeclared `Error` dominates all.
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
    // Overlapping interior grids so the product is always a Function; the right side stays Error
    // (compose(Error, Error) = Error), and only the left tail varies.
    let a = Distribution::Function(make_function(0.0, 10.0, 101, 5.0, 2.0).with_left_extrap(a_left).unwrap());
    let b = Distribution::Function(make_function(2.0, 12.0, 101, 7.0, 2.0).with_left_extrap(b_left).unwrap());

    let ab = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(fab) = &ab else {
      panic!("Expected Function, got {ab:?}")
    };
    assert_eq!(expected_left, fab.left_extrap());
    assert_eq!(BoundaryBehavior::Error, fab.right_extrap());

    // Composition is symmetric: swapping operands gives the same result tail.
    let ba = distribution_multiplication(&b, &a).unwrap();
    let Distribution::Function(fba) = &ba else {
      panic!("Expected Function, got {ba:?}")
    };
    assert_eq!(expected_left, fba.left_extrap());
  }

  /// Two soft `Linear` tails compose to a soft `Linear` tail whose neg-log slopes add
  /// (multiplication is addition in neg-log space).
  #[test]
  fn test_multiply_function_function_composes_soft_tails() {
    let a = Distribution::Function(
      make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.002 }))
        .unwrap(),
    );
    let b = Distribution::Function(
      make_function(2.0, 12.0, 101, 7.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.003 }))
        .unwrap(),
    );

    let ab = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(fab) = &ab else {
      panic!("Expected Function, got {ab:?}")
    };
    let BoundaryBehavior::Linear(law) = fab.left_extrap() else {
      panic!("Expected a composed Linear left tail, got {:?}", fab.left_extrap())
    };
    // Oracle: SoftTailLaw::compose_multiply adds the slopes: -0.002 + -0.003 = -0.005.
    assert_ulps_eq!(-0.005, law.slope, max_ulps = 4);
    assert_eq!(BoundaryBehavior::Error, fab.right_extrap());
  }

  /// Two `HardApproach` tails cannot be multiplied. At distinct boundaries the tighter bound
  /// dominates while the other operand contributes a smooth factor there, which no single-exponent
  /// hard-approach law can represent. No production path multiplies two hard-approach tails, so this
  /// is a loud internal error rather than a silent lossy composition.
  #[test]
  fn test_multiply_function_function_two_hard_approach_is_unsupported() {
    let law_a = HardApproachLaw { t_hard: 0.0, b: 1.0 };
    let law_b = HardApproachLaw { t_hard: 0.0, b: 2.0 };

    let a = Distribution::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law_a))
        .unwrap(),
    );
    let b = Distribution::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law_b))
        .unwrap(),
    );

    assert_error!(
      distribution_multiplication(&a, &b),
      "Cannot multiply two HardApproach tails: their product is not representable by a single-parameter hard-approach law, and this composition is unreachable in the inference pipeline. This is an internal error. Please report it to developers."
    );
  }

  /// A nullary `Hard` operand absorbs the other side's approach law. It declares zero density in the
  /// sub-grid gap `[t_hard, t_first)`, so the product vanishes there and the present law must not
  /// survive: the result is a nullary `Hard`, not `HardApproach`.
  #[test]
  fn test_multiply_function_function_nullary_hard_absorbs_approach_law() {
    let law = HardApproachLaw { t_hard: 0.0, b: 1.5 };

    let a = Distribution::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::HardApproach(law))
        .unwrap(),
    );
    let b = Distribution::Function(
      make_function(1.0, 10.0, 91, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Hard)
        .unwrap(),
    );

    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    assert_eq!(BoundaryBehavior::Hard, f.left_extrap());
  }

  /// normalize() preserves approach laws alongside tail policies. Under `Plain` policy normalize
  /// rescales every ordinate by `1/max`, which scales the edge-relative shape parameter by that
  /// factor so that evaluating the law on the rescaled ordinates matches rescaling the original
  /// law's output.
  #[test]
  fn test_multiply_normalize_preserves_approach_law() {
    // Build a distribution with a known max > 1 so normalization visibly rescales
    let f = DistributionFunction::<f64, Plain>::from_range_values((1.0, 5.0), array![20.0, 40.0, 100.0, 40.0, 20.0])
      .unwrap()
      .with_left_extrap(BoundaryBehavior::HardApproach(HardApproachLaw { t_hard: 0.0, b: 1.0 }))
      .unwrap();

    let normalized = Distribution::Function(f).normalize();
    let Distribution::Function(f) = &normalized else {
      panic!("Expected Function")
    };
    let preserved = f
      .left_extrap()
      .approach_law()
      .expect("approach law should survive normalization");
    // max was 100.0 -> scale factor 1/100 -> the exponent scales by 1/100
    assert_abs_diff_eq!(0.0, preserved.t_hard, epsilon = 1e-14);
    assert_abs_diff_eq!(0.01, preserved.b, epsilon = 1e-14);
  }

  /// An empty product is legitimate only when the operands' hard domains are genuinely disjoint.
  /// A soft boundary is unbounded on that side, so it bridges a gap and the domains overlap; only
  /// two hard bounds facing each other can separate the domains. Multiplication returns `Empty`
  /// exactly in the disjoint case and reports an internal error otherwise (numerical collapse).
  ///
  /// `L` = Linear (soft), `H` = Hard, `E` = Error. Each operand is (lo, hi, left_tail, right_tail).
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

  /// A multiplication that collapses to empty while the operands' hard domains overlap is a
  /// numerical failure, not a structural one, so it raises an internal error instead of silently
  /// returning `Empty`. Here the point sits inside the function's domain but lands on a zero
  /// ordinate, so the product underflows to an undefined amplitude (`Plain::is_defined` rejects it).
  #[test]
  fn test_multiply_point_on_function_zero_raises_internal_error() {
    let point = Distribution::point(1.0, 5.0);
    let func = Distribution::function(array![0.0, 1.0, 2.0], array![1.0, 0.0, 1.0]).unwrap();
    let error = distribution_multiplication(&point, &func).unwrap_err().to_string();
    assert!(error.contains("hard domains overlap"), "unexpected error: {error}");
  }

  /// Two point masses at distinct times have genuinely disjoint hard domains, so their product is a
  /// legitimate `Empty`.
  #[test]
  fn test_multiply_disjoint_points_return_empty() {
    let a = Distribution::point(1.0, 2.0);
    let b = Distribution::point(5.0, 3.0);
    assert_eq!(Distribution::Empty, distribution_multiplication(&a, &b).unwrap());
  }

  /// The single empty-result guard, shared by multiplication, convolution, and division, permits
  /// `Empty` only when the operands' hard domains are genuinely disjoint (or an operand is already
  /// empty). Overlapping hard domains mean the empty result is a numerical collapse, which raises an
  /// internal error for every operation. This is the see-red half of the invariant.
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

  /// Genuinely disjoint hard domains yield a legitimate `Empty` for every operation.
  #[rstest]
  #[case::multiplication("multiplication")]
  #[case::convolution("convolution")]
  #[case::division("division")]
  fn test_guarded_empty_result_disjoint_returns_empty(#[case] operation: &str) {
    let a = Some(((0.0, 8.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    let b = Some(((10.0, 20.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    assert_eq!(
      Distribution::Empty,
      guarded_empty_result::<Plain>(operation, a, b).unwrap()
    );
  }

  /// An operand with empty support (`None`) is disjoint from anything, so the result is a legitimate
  /// `Empty` even when the other operand covers an overlapping range.
  #[test]
  fn test_guarded_empty_result_empty_operand_returns_empty() {
    let present = Some(((0.0, 10.0), (BoundaryBehavior::Error, BoundaryBehavior::Error)));
    assert_eq!(
      Distribution::Empty,
      guarded_empty_result::<Plain>("multiplication", None, present).unwrap()
    );
    assert_eq!(
      Distribution::Empty,
      guarded_empty_result::<Plain>("multiplication", present, None).unwrap()
    );
  }
}
