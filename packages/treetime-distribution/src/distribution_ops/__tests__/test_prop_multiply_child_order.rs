//! Child-order invariance of the backward-pass message combination.
//!
//! The timetree backward pass combines a node's child parent-time messages by folding
//! `distribution_multiplication` with `normalize()` between steps
//! (`packages/treetime/src/timetree/inference/backward_pass.rs`). Each message carries a soft
//! `Linear` left tail (the parent may be arbitrarily older) and a `Hard` right tail (the
//! child's sampling date bounds the parent's age)
//! (`kb/decisions/distribution-tails-and-arithmetic.md`).
//!
//! The combined result must not depend on the order children are processed. Multiplication is
//! commutative, and both the support intersection and the composed tails are built from
//! order-independent `max`/`min` over the operands' bounds, so across every permutation of the
//! children the fold yields a non-empty distribution over the *same* support with the *same*
//! `Linear`/`Hard` tails. These invariants guard the exact failure modes the fix targeted:
//! tail-metadata loss (which would collapse some orderings to `Empty` while others survive) and
//! asymmetric support resolution (which would place the result on different intervals).
//!
//! Scope of the invariant: the *support and tails* are order-invariant; the fine-grained grid
//! resolution and hence the argmax `likely_time()` are not. Sequential multiplication resamples
//! at each binary step and is not associative -- the point count `round(width / dx) + 1` drifts
//! with the order-dependent intermediate widths, so peaks agree only to grid resolution
//! (`kb/decisions/distribution-tails-and-arithmetic.md`, "Accepted limitations"). This suite
//! therefore asserts the exact structural invariants and does not assert `likely_time()`
//! equality across orderings.

#[cfg(test)]
mod tests {
  use crate::DistributionFunction;
  use crate::DistributionPlain as Distribution;
  use crate::distribution_ops::multiply::distribution_multiplication;
  use crate::policy::Plain;
  use itertools::Itertools;
  use ndarray::Array1;
  use proptest::prelude::*;
  use treetime_grid::{BoundaryBehavior, SoftTailLaw};
  use treetime_utils::{pretty_assert_ulps_eq, prop_assert_ulps_eq};

  /// A gentle soft `Linear` left tail (the backward-pass parent-older tail). `Linear` is the only
  /// soft `BoundaryBehavior`, so it is what bridges disjoint grids; the slope is tiny and negative so
  /// the extrapolated Plain ordinates stay positive across the wide gaps these fixtures span.
  const SOFT: BoundaryBehavior = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.001 });

  proptest! {
    /// Every permutation of the children folds to a non-empty `Function` over the same support.
    ///
    /// A lost soft left tail would let a disjoint child collapse the product to `Empty` for
    /// some orderings but not others; an asymmetric intersection would shift the support bounds.
    #[test]
    fn test_prop_multiply_child_order_preserves_support(
      params in prop::collection::vec(child_message_strategy(), 2..=5),
    ) {
      let children: Vec<Distribution> = params.iter().map(|&p| make_child_message(p)).collect();
      let children: Vec<&Distribution> = children.iter().collect();
      let n = children.len();

      let (base_min, base_max) = fold_support(&children)?;
      for perm in children.iter().copied().permutations(n) {
        let (perm_min, perm_max) = fold_support(&perm)?;
        // Support bounds come from order-independent max/min over the operands' grid bounds, so
        // they are the same for every ordering (max_ulps guards resampling round-off only).
        prop_assert_ulps_eq!(base_min, perm_min, max_ulps = 4);
        prop_assert_ulps_eq!(base_max, perm_max, max_ulps = 4);
      }
    }

    /// Every permutation of the children folds to a result carrying `Linear` left / `Hard` right
    /// tails, so a later disjoint child remains reachable regardless of processing order.
    #[test]
    fn test_prop_multiply_child_order_preserves_tails(
      params in prop::collection::vec(child_message_strategy(), 2..=5),
    ) {
      let children: Vec<Distribution> = params.iter().map(|&p| make_child_message(p)).collect();
      let children: Vec<&Distribution> = children.iter().collect();
      let n = children.len();

      for perm in children.iter().copied().permutations(n) {
        let result = fold_children(&perm);
        let Distribution::Function(f) = &result else {
          return Err(TestCaseError::fail(format!("fold collapsed to {result:?}")));
        };
        prop_assert!(matches!(f.left_extrap(), BoundaryBehavior::Linear(_)));
        prop_assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
      }
    }
  }

  /// A disjoint child reachable only via the soft `Linear` left tail keeps the fold non-empty and on
  /// the same support no matter where it appears in the processing order.
  ///
  /// Two recent messages narrow the accumulator to a late interval; a much older message overlaps
  /// only through its soft `Linear` left tail. Before tail-preserving normalization, folding the old
  /// message last collapsed the product to `Empty`; the support is now identical for every order.
  #[test]
  fn test_multiply_child_order_disjoint_child_invariant_support() {
    let recent_a = make_child_message((2025.0, 2.0, 0.5, 41));
    let recent_b = make_child_message((2024.5, 1.0, 0.3, 21));
    let old = make_child_message((1985.0, 30.0, 5.0, 291));

    let children = [&recent_a, &recent_b, &old];
    let (base_min, base_max) = fold_support(&children).unwrap();

    for perm in children.iter().copied().permutations(3) {
      let (perm_min, perm_max) = fold_support(&perm).unwrap();
      pretty_assert_ulps_eq!(base_min, perm_min, max_ulps = 4);
      pretty_assert_ulps_eq!(base_max, perm_max, max_ulps = 4);
    }
  }

  /// Fold the children with multiply + normalize (the backward-pass combination) and return the
  /// result's support bounds, failing the property if any ordering collapses to a non-`Function`.
  fn fold_support(children: &[&Distribution]) -> Result<(f64, f64), TestCaseError> {
    let result = fold_children(children);
    let Distribution::Function(f) = &result else {
      return Err(TestCaseError::fail(format!("fold collapsed to {result:?}")));
    };
    Ok((f.x_min(), f.x_max()))
  }

  /// Combine child messages exactly as `propagate_distributions_backward_slot` does: multiply the
  /// running accumulator by each message and normalize between steps.
  fn fold_children(children: &[&Distribution]) -> Distribution {
    let mut accum = children[0].clone();
    for child in &children[1..] {
      accum = distribution_multiplication(&accum, child).unwrap().normalize();
    }
    accum
  }

  /// A backward parent message: a Gaussian bump on a finite grid with a soft `Linear` left tail
  /// (parent may be older) and a `Hard` right tail (child date bounds the parent's age).
  fn make_child_message((center, width, sigma, n_points): (f64, f64, f64, usize)) -> Distribution {
    let x_min = center - width / 2.0;
    let dx = width / (n_points - 1) as f64;
    let y = Array1::from_shape_fn(n_points, |i| {
      let x = x_min + dx * i as f64;
      (-0.5 * ((x - center) / sigma).powi(2)).exp()
    });
    let f = DistributionFunction::<f64, Plain>::from_start_dx_values(x_min, dx, y)
      .unwrap()
      .with_left_extrap(SOFT)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();
    Distribution::Function(f)
  }

  /// `(center, width, sigma, n_points)` for a child message. The wide center range lets supports
  /// be disjoint so tail-driven overlap is exercised, not only genuine grid overlap.
  fn child_message_strategy() -> impl Strategy<Value = (f64, f64, f64, usize)> {
    (1985.0_f64..2015.0, 2.0_f64..30.0, 0.5_f64..8.0, 20_usize..200)
  }
}
