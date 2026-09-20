#[cfg(test)]
mod tests {
  use crate::piecewise_linear_fn::PiecewiseLinearFn;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_piecewise_linear_eval_interpolation() {
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

    pretty_assert_ulps_eq!(pl.eval(-5.0), 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(-100.0), 0.0, max_ulps = 4);

    pretty_assert_ulps_eq!(pl.eval(15.0), 100.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(1000.0), 100.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_linear_eval_multiple_segments() {
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0, 20.0, 30.0], array![0.0, 50.0, 100.0, 0.0]);

    pretty_assert_ulps_eq!(pl.eval(5.0), 25.0, max_ulps = 4);

    pretty_assert_ulps_eq!(pl.eval(15.0), 75.0, max_ulps = 4);

    pretty_assert_ulps_eq!(pl.eval(25.0), 50.0, max_ulps = 4);

    pretty_assert_ulps_eq!(pl.eval(10.0), 50.0, max_ulps = 4);
    pretty_assert_ulps_eq!(pl.eval(20.0), 100.0, max_ulps = 4);
  }

  #[test]
  fn test_piecewise_linear_eval_many() {
    let pl = PiecewiseLinearFn::new(array![0.0, 10.0], array![0.0, 100.0]);
    let ts = array![-5.0, 0.0, 5.0, 10.0, 15.0];
    let result = pl.eval_many(&ts);

    pretty_assert_ulps_eq!(result[0], 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[1], 0.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[2], 50.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[3], 100.0, max_ulps = 4);
    pretty_assert_ulps_eq!(result[4], 100.0, max_ulps = 4);
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
