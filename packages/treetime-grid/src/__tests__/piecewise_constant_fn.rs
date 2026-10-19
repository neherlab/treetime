#[cfg(test)]
mod tests {
  use crate::piecewise_constant_fn::PiecewiseConstantFn;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_piecewise_constant_eval() {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_ulps_eq!(pc.eval(-1.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(0.5), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(1.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(3.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(5.0), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(7.5), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(10.0), 3.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval(100.0), 3.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_constant_eval_left() {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    pretty_assert_ulps_eq!(pc.eval_left(-1.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(0.5), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(1.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(3.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(5.0), 1.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(7.5), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(10.0), 2.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pc.eval_left(100.0), 3.0, max_ulps = 4);
  }
}
