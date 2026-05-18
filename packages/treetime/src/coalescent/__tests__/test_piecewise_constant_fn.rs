#[cfg(test)]
mod tests {
  use crate::coalescent::piecewise_constant_fn::PiecewiseConstantFn;
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

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
