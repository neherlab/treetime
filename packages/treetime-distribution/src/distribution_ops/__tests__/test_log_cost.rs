#[cfg(test)]
mod tests {
  use crate::DistributionFunction;
  use crate::DistributionNegLog;
  use crate::distribution_ops::log_cost::distribution_add_neg_log_weight;
  use approx::assert_abs_diff_eq;
  use ndarray::{Array1, array};
  use treetime_grid::{BoundaryBehavior, DEFAULT_TAIL_FIT_POINTS, Side, SoftTailLaw};
  use treetime_utils::assert_error;

  /// A constant weight is a uniform shift in log space, which normalization removes: the result
  /// equals the already-normalized input.
  #[test]
  fn test_log_cost_add_neg_log_weight_constant_is_uniform_shift() {
    let neglog = DistributionNegLog::function(array![0.0, 1.0, 2.0], array![4.0, 0.0, 3.0]).unwrap();
    let actual = distribution_add_neg_log_weight(&neglog, |_| Ok(10.0)).unwrap();
    assert_abs_diff_eq!(array![4.0, 0.0, 3.0], actual.y(), epsilon = 1e-15);
  }

  /// Empty passes through unchanged.
  #[test]
  fn test_log_cost_add_neg_log_weight_empty_passthrough() {
    let actual = distribution_add_neg_log_weight(&DistributionNegLog::Empty, |_| Ok(1.0)).unwrap();
    assert_eq!(DistributionNegLog::Empty, actual);
  }

  /// A Formula has no grid and is rejected.
  #[test]
  fn test_log_cost_add_neg_log_weight_rejects_formula() {
    let formula = DistributionNegLog::Formula(crate::DistributionFormula::new(|_| Ok(1.0), 0.0, 10.0));
    assert_error!(
      distribution_add_neg_log_weight(&formula, |_| Ok(1.0)),
      "distribution_add_neg_log_weight requires a concrete Point, Range, or Function distribution"
    );
  }

  /// A varying weight keeps the input's per-side tail *class* and re-fits the soft slope it changed.
  ///
  /// The input is `y = -2t` with a fitted soft `Linear` left tail (slope `-2`) and a `Hard` right
  /// edge -- the exact Linear-left / Hard-right policy the coalescent backward pass produces. Adding
  /// `w(t) = 0.5 t` makes the combined ordinate `-1.5 t`, so the re-fit left slope is `-1.5`
  /// (analytical: least-squares slope of an exactly linear ordinate) and the right edge stays `Hard`.
  /// Before the fix the rebuild reset both sides to `Error`, disabling the following mass re-window.
  #[test]
  fn test_log_cost_add_neg_log_weight_preserves_soft_left_hard_right() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let f = DistributionFunction::from_arrays(&t, y).unwrap();
    let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS).unwrap();
    let input = f
      .with_left_extrap(BoundaryBehavior::Linear(left))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_add_neg_log_weight(&DistributionNegLog::Function(input), |ti: f64| Ok(0.5 * ti)).unwrap();

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

  /// A constant weight is a pure shift, so the re-fit leaves the soft slope unchanged.
  ///
  /// Same `y = -2t` input (left slope `-2`); a constant weight adds no `t`-dependence, so the
  /// combined ordinate is still `-2t` and the re-fit recovers slope `-2`. This is the boundary case
  /// where the general re-fit reduces to `normalize`'s shift-invariant carry.
  #[test]
  fn test_log_cost_add_neg_log_weight_constant_weight_keeps_soft_slope() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let f = DistributionFunction::from_arrays(&t, y).unwrap();
    let left = SoftTailLaw::fit(f.grid_fn(), Side::Left, DEFAULT_TAIL_FIT_POINTS).unwrap();
    let input = f
      .with_left_extrap(BoundaryBehavior::Linear(left))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_add_neg_log_weight(&DistributionNegLog::Function(input), |_| Ok(10.0)).unwrap();

    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result");
    };
    let BoundaryBehavior::Linear(law) = rf.left_extrap() else {
      panic!("left tail must stay soft Linear, got {:?}", rf.left_extrap());
    };
    assert_abs_diff_eq!(-2.0, law.slope, epsilon = 1e-12);
    assert_eq!(BoundaryBehavior::Hard, rf.right_extrap());
  }

  /// Hard sides carry through unchanged: `0 * exp(-w) = 0` beyond the edge regardless of the weight.
  #[test]
  fn test_log_cost_add_neg_log_weight_carries_hard_both_sides() {
    let t = Array1::linspace(0.0, 4.0, 21);
    let y = t.mapv(|ti| -2.0 * ti);
    let input = DistributionFunction::from_arrays(&t, y)
      .unwrap()
      .with_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_add_neg_log_weight(&DistributionNegLog::Function(input), |ti: f64| Ok(0.5 * ti)).unwrap();

    let DistributionNegLog::Function(rf) = actual else {
      panic!("expected a Function result");
    };
    assert_eq!(BoundaryBehavior::Hard, rf.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, rf.right_extrap());
  }
}
