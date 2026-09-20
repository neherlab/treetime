#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

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

  const SOFT: BoundaryBehavior = BoundaryBehavior::Linear(SoftTailLaw { slope: -0.001 });

  proptest! {
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
        prop_assert_ulps_eq!(base_min, perm_min, max_ulps = 4);
        prop_assert_ulps_eq!(base_max, perm_max, max_ulps = 4);
      }
    }

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

  fn fold_support(children: &[&Distribution]) -> Result<(f64, f64), TestCaseError> {
    let result = fold_children(children);
    let Distribution::Function(f) = &result else {
      return Err(TestCaseError::fail(format!("fold collapsed to {result:?}")));
    };
    Ok((f.x_min(), f.x_max()))
  }

  fn fold_children(children: &[&Distribution]) -> Distribution {
    let mut accum = children[0].clone();
    for child in &children[1..] {
      accum = distribution_multiplication(&accum, child).unwrap().normalize();
    }
    accum
  }

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

  fn child_message_strategy() -> impl Strategy<Value = (f64, f64, f64, usize)> {
    (1985.0_f64..2015.0, 2.0_f64..30.0, 0.5_f64..8.0, 20_usize..200)
  }
}
