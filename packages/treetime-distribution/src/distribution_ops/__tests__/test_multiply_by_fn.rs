#[cfg(test)]
mod tests {
  use crate::__tests__::aliases::DistributionNegLog;
  use crate::DistributionFunction;
  use crate::distribution_ops::multiply_by_fn::distribution_multiply_by_fn;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array1, array};
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};
  use treetime_utils::assert_error;

  #[test]
  fn test_multiply_by_fn_constant_is_uniform_shift() {
    let neglog = DistributionNegLog::function(array![0.0, 1.0, 2.0], array![4.0, 0.0, 3.0]).unwrap();
    let actual = distribution_multiply_by_fn(&neglog, |_| Ok(10.0)).unwrap();
    assert_abs_diff_eq!(array![4.0, 0.0, 3.0], actual.y().unwrap(), epsilon = 1e-15);
  }

  #[test]
  fn test_multiply_by_fn_empty_passthrough() {
    let actual = distribution_multiply_by_fn(&DistributionNegLog::Empty, |_| Ok(1.0)).unwrap();
    assert_eq!(DistributionNegLog::Empty, actual);
  }

  #[test]
  fn test_multiply_by_fn_rejects_formula() {
    let formula = DistributionNegLog::Formula(crate::DistributionFormula::new(|_| Ok(1.0), 0.0, 10.0));
    assert_error!(
      distribution_multiply_by_fn(&formula, |_| Ok(1.0)),
      "distribution_multiply_by_fn requires a concrete Point, Range, or Function distribution"
    );
  }

  #[test]
  fn test_multiply_by_fn_preserves_soft_left_hard_right() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let f = DistributionFunction::from_arrays(&t, y).unwrap();
    let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS).unwrap();
    let input = f
      .with_left_extrap(BoundaryBehavior::Linear(left))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_multiply_by_fn(&DistributionNegLog::Function(input), |ti: f64| Ok(0.5 * ti)).unwrap();

    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result");
    };
    let BoundaryBehavior::Linear(law) = rf.left_extrap() else {
      panic!(
        "left tail must stay soft Linear, not reset to Error, got {:?}",
        rf.left_extrap()
      );
    };
    assert_abs_diff_eq!(-1.5, law.slope, epsilon = 1e-12);
    assert_eq!(BoundaryBehavior::Hard, rf.right_extrap());
  }

  #[test]
  fn test_multiply_by_fn_constant_weight_keeps_soft_slope() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let f = DistributionFunction::from_arrays(&t, y).unwrap();
    let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS).unwrap();
    let input = f
      .with_left_extrap(BoundaryBehavior::Linear(left))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_multiply_by_fn(&DistributionNegLog::Function(input), |_| Ok(10.0)).unwrap();

    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result");
    };
    let BoundaryBehavior::Linear(law) = rf.left_extrap() else {
      panic!("left tail must stay soft Linear, got {:?}", rf.left_extrap());
    };
    assert_abs_diff_eq!(-2.0, law.slope, epsilon = 1e-12);
    assert_eq!(BoundaryBehavior::Hard, rf.right_extrap());
  }

  #[test]
  fn test_multiply_by_fn_carries_hard_both_sides() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let input = DistributionFunction::from_arrays(&t, y)
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Hard)
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_multiply_by_fn(&DistributionNegLog::Function(input), |ti: f64| Ok(0.5 * ti)).unwrap();

    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result");
    };
    assert_eq!(BoundaryBehavior::Hard, rf.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, rf.right_extrap());
  }
}
