#[cfg(test)]
mod tests {
  use crate::piecewise_constant_fn::PiecewiseConstantFn;
  use ndarray::array;
  use rstest::rstest;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::below_first_breakpoint(  -1.0, 0.0)]
  #[case::at_first_breakpoint(      1.0, 1.0)]
  #[case::inside_an_interval(       3.0, 1.0)]
  #[case::at_last_breakpoint(      10.0, 3.0)]
  #[case::above_last_breakpoint(  100.0, 3.0)]
  #[trace]
  fn test_piecewise_constant_eval_is_right_continuous(#[case] t: f64, #[case] expected: f64) {
    let pc = PiecewiseConstantFn::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);

    let actual = pc.eval(t);

    pretty_assert_ulps_eq!(expected, actual, max_ulps = 4);
  }
}
