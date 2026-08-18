#[cfg(test)]
mod tests {
  use crate::DistributionPlain as Distribution;
  use crate::distribution_core::function::DistributionFunction;
  use crate::distribution_ops::subtract::distribution_subtraction;
  use approx::assert_ulps_eq;
  use ndarray::array;
  use treetime_grid::{BoundaryBehavior, SoftTailLaw};
  use treetime_utils::assert_error;

  /// A pointwise difference on matching grids subtracts ordinates and keeps `Hard` on both sides when
  /// both operands vanish beyond the edge (`0 - 0 = 0`).
  #[test]
  fn test_subtract_function_hard_both_sides_stays_hard() {
    let a = DistributionFunction::from_start_dx_values(0.0, 1.0, array![5.0, 4.0, 3.0])
      .unwrap()
      .with_extrap(BoundaryBehavior::Hard)
      .unwrap();
    let b = DistributionFunction::from_start_dx_values(0.0, 1.0, array![1.0, 1.0, 1.0])
      .unwrap()
      .with_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_subtraction(&Distribution::Function(a), &Distribution::Function(b)).unwrap();

    let Distribution::Function(f) = actual else {
      panic!("expected a Function difference");
    };
    assert_ulps_eq!(f.y(), &array![4.0, 3.0, 2.0], max_ulps = 4);
    assert_eq!(BoundaryBehavior::Hard, f.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
  }

  /// A soft or `Error` operand on a side leaves the difference without a representable tail law, so
  /// that side is `Error`; a side where both operands are zero-beyond stays `Hard`.
  #[test]
  fn test_subtract_function_soft_side_becomes_error() {
    let a = DistributionFunction::from_start_dx_values(0.0, 1.0, array![5.0, 4.0, 3.0])
      .unwrap()
      .with_left_extrap(BoundaryBehavior::Linear(SoftTailLaw { slope: -0.3 }))
      .unwrap()
      .with_right_extrap(BoundaryBehavior::Hard)
      .unwrap();
    let b = DistributionFunction::from_start_dx_values(0.0, 1.0, array![1.0, 1.0, 1.0])
      .unwrap()
      .with_extrap(BoundaryBehavior::Hard)
      .unwrap();

    let actual = distribution_subtraction(&Distribution::Function(a), &Distribution::Function(b)).unwrap();

    let Distribution::Function(f) = actual else {
      panic!("expected a Function difference");
    };
    assert_eq!(BoundaryBehavior::Error, f.left_extrap());
    assert_eq!(BoundaryBehavior::Hard, f.right_extrap());
  }

  /// Subtraction requires identical grids.
  #[test]
  fn test_subtract_function_mismatched_grids_errors() {
    let a = Distribution::function(array![0.0, 1.0, 2.0], array![5.0, 4.0, 3.0]).unwrap();
    let b = Distribution::function(array![0.0, 1.0, 2.0, 3.0], array![1.0, 1.0, 1.0, 1.0]).unwrap();

    assert_error!(
      distribution_subtraction(&a, &b),
      "Cannot subtract distributions with different grids"
    );
  }

  /// Non-`Function` operands are rejected.
  #[test]
  fn test_subtract_non_function_errors() {
    let a = Distribution::function(array![0.0, 1.0, 2.0], array![5.0, 4.0, 3.0]).unwrap();
    let point = Distribution::point(1.0, 2.0);

    assert_error!(
      distribution_subtraction(&a, &point),
      "Subtraction only supported for Function distributions with matching grids"
    );
  }
}
