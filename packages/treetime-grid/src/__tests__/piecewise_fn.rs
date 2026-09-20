#[cfg(test)]
mod tests {
  use crate::piecewise_fn::PiecewiseFnBase;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_piecewise_fn_base_accessors() {
    let base = PiecewiseFnBase::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);
    pretty_assert_ulps_eq!(&array![1.0, 5.0, 10.0], base.breakpoints(), max_ulps = 4);
    pretty_assert_ulps_eq!(&array![0.0, 1.0, 2.0, 3.0], base.values(), max_ulps = 4);
    pretty_assert_ulps_eq!(&[1.0, 5.0, 10.0][..], base.breakpoints_slice(), max_ulps = 4);
    pretty_assert_ulps_eq!(&[0.0, 1.0, 2.0, 3.0][..], base.values_slice(), max_ulps = 4);
  }

  #[test]
  fn test_piecewise_fn_base_single_breakpoint() {
    let base = PiecewiseFnBase::new(array![5.0], array![0.0, 1.0]);
    assert_eq!(1, base.breakpoints().len());
    assert_eq!(2, base.values().len());
  }

  #[test]
  fn test_piecewise_fn_base_clone() {
    let base = PiecewiseFnBase::new(array![1.0, 2.0], array![10.0, 20.0]);
    let cloned = base.clone();
    assert_eq!(cloned.breakpoints(), base.breakpoints());
    assert_eq!(cloned.values(), base.values());
  }
}
