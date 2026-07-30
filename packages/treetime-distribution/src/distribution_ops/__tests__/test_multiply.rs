#[cfg(test)]
mod tests {
  use crate::DistributionFunction;
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::formula::DistributionFormula;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::policy::Plain;
  use approx::assert_ulps_eq;
  use ndarray::{Array1, array};
  use rstest::rstest;
  use treetime_grid::BoundaryBehavior;
  use treetime_utils::pretty_assert_ulps_eq;

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

  /// C2: overlapping grids with Constant tails -- intersection unchanged (tails only extend).
  #[test]
  fn test_multiply_tail_c2_overlapping_with_tails() {
    let a = Distribution::Function(
      make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap(),
    );
    let b = Distribution::Function(
      make_function(3.0, 13.0, 101, 8.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap(),
    );
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // A's Constant left extends to min(0, 3) = 0 (no change). B's extends to min(3, 0) = 0.
    // Intersection of [0, 10] and [0, 13] = [0, 10].
    assert_ulps_eq!(f.x_min(), 0.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 10.0, max_ulps = 4);
  }

  /// C3: disjoint grids, one Constant left tail -- extends to produce non-empty result.
  #[test]
  fn test_multiply_tail_c3_disjoint_one_constant_left() {
    // A: grid [10, 20], Constant left tail
    // B: grid [0, 8], default Error tails
    let a = Distribution::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(BoundaryBehavior::Constant)
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

  /// C4: disjoint grids, both Constant left -- backward-pass pattern.
  #[test]
  fn test_multiply_tail_c4_disjoint_both_constant_left() {
    // Two backward-pass messages with disjoint finite grids but overlapping via Constant
    // left tails: leaf message [2001, 2007] and subtree message [1970, 2000].
    let leaf_msg = Distribution::Function(
      make_function(2001.0, 2007.0, 61, 2004.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let subtree_msg = Distribution::Function(
      make_function(1970.0, 2000.0, 301, 1990.0, 5.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let result = distribution_multiplication(&leaf_msg, &subtree_msg).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function (non-empty product), got {result:?}")
    };
    // Leaf extends left to min(2001, 1970) = 1970. Subtree extends left to min(1970, 2001) = 1970.
    // Right sides: both Zero, no extension.
    // Intersection of [1970, 2007] and [1970, 2000] = [1970, 2000].
    assert_ulps_eq!(f.x_min(), 1970.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 2000.0, max_ulps = 4);
    // Product should be dominated by the subtree's shape (leaf is flat constant in this region).
    let peak = result.likely_time().unwrap();
    assert!(peak > 1985.0 && peak < 1995.0, "Peak at {peak}, expected near 1990");
  }

  /// C5: disjoint grids, Zero tail -- stays Empty (zero outside support, intersection correct).
  #[test]
  fn test_multiply_tail_c5_disjoint_zero_tail() {
    let a = Distribution::Function(
      make_function(10.0, 20.0, 101, 15.0, 3.0)
        .with_left_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let b = Distribution::Function(make_function(0.0, 8.0, 101, 4.0, 2.0));
    let result = distribution_multiplication(&a, &b).unwrap();
    assert!(
      matches!(result, Distribution::Empty),
      "Zero tail should not prevent Empty"
    );
  }

  /// C6: endpoint contact with Constant tail -- extends past contact to produce interval.
  #[test]
  fn test_multiply_tail_c6_endpoint_contact_with_constant() {
    // Without tails: endpoint contact -> Point. With Constant left on A: extends past.
    let a = Distribution::Function(
      make_function(5.0, 10.0, 51, 7.5, 2.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap(),
    );
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
    let outer = Distribution::Function(
      make_function(0.0, 20.0, 201, 10.0, 5.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap(),
    );
    let inner = Distribution::Function(make_function(5.0, 15.0, 101, 10.0, 3.0));
    let result = distribution_multiplication(&outer, &inner).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // Inner is fully contained -- tails don't change the intersection.
    assert_ulps_eq!(f.x_min(), 5.0, max_ulps = 4);
    assert_ulps_eq!(f.x_max(), 15.0, max_ulps = 4);
  }

  /// C8: mixed Constant left + Zero right -- backward message pattern.
  #[test]
  fn test_multiply_tail_c8_mixed_constant_zero() {
    // Both messages have Constant left (parent could be older) and Zero right (hard upper bound).
    // Overlapping case -- verify tails don't corrupt the result.
    let a = Distribution::Function(
      make_function(0.0, 10.0, 101, 5.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let b = Distribution::Function(
      make_function(3.0, 8.0, 51, 5.5, 1.5)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let result = distribution_multiplication(&a, &b).unwrap();
    let Distribution::Function(f) = &result else {
      panic!("Expected Function")
    };
    // Both extend left to min of both = 0. Right: both Zero, no extension.
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
        .with_left_extrap(BoundaryBehavior::Constant)
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

  /// Chained multiplication with normalize-then-reapply-tails survives disjoint late child.
  ///
  /// Simulates the backward pass: three overlapping messages produce a narrow accumulated
  /// result. After normalize (which resets tails to Error), the accumulated result must have
  /// tails re-applied before multiplying with a fourth disjoint message. Without re-application,
  /// the fourth multiplication collapses to Empty.
  #[test]
  fn test_multiply_tail_chained_normalize_reapply() {
    let msg1 = Distribution::Function(
      make_function(2000.0, 2010.0, 101, 2005.0, 2.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let msg2 = Distribution::Function(
      make_function(2001.0, 2008.0, 71, 2004.0, 1.5)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let msg3 = Distribution::Function(
      make_function(2003.0, 2009.0, 61, 2006.0, 1.5)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    // Fourth message: disjoint from the narrowed accumulator but reachable via Constant left tail.
    let msg4 = Distribution::Function(
      make_function(1970.0, 1999.0, 291, 1990.0, 5.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );

    // Simulate backward pass accumulation with normalize + tail re-application.
    let mut accum = msg1;
    for msg in [&msg2, &msg3, &msg4] {
      accum = distribution_multiplication(&accum, msg)
        .unwrap()
        .normalize()
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap();
    }

    assert!(
      !matches!(accum, Distribution::Empty),
      "Chained multiplication collapsed to Empty despite Constant left tails"
    );
    assert!(
      accum.likely_time().is_some(),
      "Accumulated result must have a likely_time"
    );
  }

  /// Chained multiplication where normalize resets the accumulator's tails to Error.
  ///
  /// Two recent messages narrow the accumulator to ~(2024.5, 2025.5). After normalize resets
  /// tails to Error, a much older message at (1970, 1999) is disjoint. The old message's
  /// Constant left can't help (it extends leftward, not rightward). The accumulator's Error
  /// left can't extend leftward either. Product is Empty -- this is the scenario the backward
  /// pass re-apply fix addresses.
  #[test]
  fn test_multiply_tail_chained_normalize_without_reapply_collapses() {
    let msg_recent_1 = Distribution::Function(
      make_function(2024.0, 2026.0, 21, 2025.0, 0.5)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let msg_recent_2 = Distribution::Function(
      make_function(2024.5, 2025.5, 11, 2025.0, 0.3)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );
    let msg_old = Distribution::Function(
      make_function(1970.0, 1999.0, 291, 1990.0, 5.0)
        .with_left_extrap(BoundaryBehavior::Constant)
        .unwrap()
        .with_right_extrap(BoundaryBehavior::Zero)
        .unwrap(),
    );

    // First two messages overlap and narrow the accumulator.
    let product = distribution_multiplication(&msg_recent_1, &msg_recent_2).unwrap();
    assert!(!matches!(product, Distribution::Empty));

    // Normalize resets tails to Error.
    let normalized = product.normalize();

    // Old message's Constant left extends to min(1970, ~2024.5) = 1970 (already there).
    // Accumulator's Error left stays at ~2024.5. Disjoint.
    let final_result = distribution_multiplication(&normalized, &msg_old).unwrap();
    assert!(
      matches!(final_result, Distribution::Empty),
      "Without tail re-application, the accumulator can't extend leftward to overlap the old message"
    );
  }
}
