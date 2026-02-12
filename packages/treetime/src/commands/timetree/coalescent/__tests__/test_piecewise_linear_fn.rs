#[cfg(test)]
mod tests {
  use crate::commands::timetree::coalescent::piecewise_linear_fn::PiecewiseLinearFn;
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_piecewise_linear_eval_interpolation() {
    // Breakpoints: [0.0, 10.0]
    // Values: [0.0, 100.0]
    // Linear from 0 to 100
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);

    pretty_assert_ulps_eq!(pl.eval(0.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(5.0), 50.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(10.0), 100.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(2.5), 25.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(7.5), 75.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_linear_eval_extrapolation() {
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);

    // Before first breakpoint: constant at first value
    pretty_assert_ulps_eq!(pl.eval(-5.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(-100.0), 0.0, max_ulps = 4);

    // After last breakpoint: constant at last value
    pretty_assert_ulps_eq!(pl.eval(15.0), 100.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(1000.0), 100.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_linear_eval_multiple_segments() {
    // Three segments: 0->10: 0->50, 10->20: 50->100, 20->30: 100->0
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0, 20.0, 30.0], array![0.0, 50.0, 100.0, 0.0]);

    // First segment
    pretty_assert_ulps_eq!(pl.eval(5.0), 25.0, max_ulps = 4);

    // Second segment
    pretty_assert_ulps_eq!(pl.eval(15.0), 75.0, max_ulps = 4);

    // Third segment (decreasing)
    pretty_assert_ulps_eq!(pl.eval(25.0), 50.0, max_ulps = 4);

    // Exact breakpoints
    pretty_assert_ulps_eq!(pl.eval(10.0), 50.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(20.0), 100.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_linear_eval_many() {
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);
    let ts = array![-5.0, 0.0, 5.0, 10.0, 15.0];
    let result = pl.eval_many(&ts);

    pretty_assert_ulps_eq!(result[0], 0.0, max_ulps = 4); // extrapolation before
    pretty_assert_ulps_eq!(result[1], 0.0, max_ulps = 4); // at first breakpoint
    pretty_assert_ulps_eq!(result[2], 50.0, max_ulps = 4); // interpolation
    pretty_assert_ulps_eq!(result[3], 100.0, max_ulps = 4); // at last breakpoint
    pretty_assert_ulps_eq!(result[4], 100.0, max_ulps = 4); // extrapolation after
  }

  #[test]
  fn test_piecewise_linear_accessors() {
    let breakpoints = array![1.0, 5.0, 10.0];
    let values = array![10.0, 20.0, 30.0];
    let pl = PiecewiseLinearFn::new(breakpoints.clone(), values.clone());

    assert_eq!(pl.breakpoints(), &breakpoints);
    assert_eq!(pl.values(), &values);
  }
}
