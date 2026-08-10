#[cfg(test)]
mod tests {
  use crate::piecewise_constant_fn::PiecewiseConstantFn;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_piecewise_constant_eval() {
    // Breakpoints: [1.0, 5.0, 10.0]
    // Values: [0.0, 1.0, 2.0, 3.0]
    // t < 1.0 -> 0.0
    // 1.0 <= t < 5.0 -> 1.0
    // 5.0 <= t < 10.0 -> 2.0
    // t >= 10.0 -> 3.0
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_ulps_eq!(pc.eval(-1.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(0.5), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(1.0), 1.0, max_ulps = 4); // at breakpoint: value after
    pretty_assert_ulps_eq!(pc.eval(3.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(5.0), 2.0, max_ulps = 4); // at breakpoint: value after
    pretty_assert_ulps_eq!(pc.eval(7.5), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(10.0), 3.0, max_ulps = 4); // at breakpoint: value after
    pretty_assert_ulps_eq!(pc.eval(100.0), 3.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_constant_eval_left() {
    // Same function as test_piecewise_constant_eval, but eval_left
    // returns the pre-event (left-limit) value at breakpoints.
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_ulps_eq!(pc.eval_left(-1.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(0.5), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(1.0), 0.0, max_ulps = 4); // at breakpoint: value BEFORE
    pretty_assert_ulps_eq!(pc.eval_left(3.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(5.0), 1.0, max_ulps = 4); // at breakpoint: value BEFORE
    pretty_assert_ulps_eq!(pc.eval_left(7.5), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(10.0), 2.0, max_ulps = 4); // at breakpoint: value BEFORE
    pretty_assert_ulps_eq!(pc.eval_left(100.0), 3.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_constant_zip_map_uses_breakpoint_union() {
    let left = PiecewiseConstantFn::new(array![1.0, 3.0], array![2.0, 4.0, 8.0]);
    let right = PiecewiseConstantFn::new(array![2.0, 3.0], array![10.0, 20.0, 40.0]);

    let actual = left.zip_map(&right, |left, right| left * right);

    pretty_assert_ulps_eq!(actual.breakpoints(), &array![1.0, 2.0, 3.0], max_ulps = 4);
    pretty_assert_ulps_eq!(actual.values(), &array![20.0, 40.0, 80.0, 320.0], max_ulps = 4);
  }

  #[test]
  fn test_piecewise_constant_zip_map_supports_constant_input() {
    let constant = PiecewiseConstantFn::new(array![], array![2.0]);
    let stepped = PiecewiseConstantFn::new(array![1.0], array![3.0, 5.0]);

    let actual = constant.zip_map(&stepped, |left, right| left + right);

    pretty_assert_ulps_eq!(actual.breakpoints(), &array![1.0], max_ulps = 4);
    pretty_assert_ulps_eq!(actual.values(), &array![5.0, 7.0], max_ulps = 4);
  }

  #[test]
  fn test_piecewise_constant_eval_many() {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0], array![0.0, 1.0, 2.0]);
    let ts = array![0.0, 1.0, 3.0, 5.0, 10.0];
    let result = pc.eval_many(&ts);
    pretty_assert_ulps_eq!(result[0], 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[1], 1.0, max_ulps = 4); // at breakpoint: value after
    pretty_assert_ulps_eq!(result[2], 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[3], 2.0, max_ulps = 4); // at breakpoint: value after
    pretty_assert_ulps_eq!(result[4], 2.0, max_ulps = 4);
  }
}
