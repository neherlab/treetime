#[cfg(test)]
mod tests {
  use crate::piecewise_fn::PiecewiseFnBase;
  use ndarray::array;

  #[test]
  fn test_piecewise_fn_base_accessors() {
    let base = PiecewiseFnBase::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);
    assert_eq!(base.breakpoints(), &array![1.0, 5.0, 10.0]);
    assert_eq!(base.values(), &array![0.0, 1.0, 2.0, 3.0]);
    assert_eq!(&[1.0, 5.0, 10.0], base.breakpoints_slice());
    assert_eq!(&[0.0, 1.0, 2.0, 3.0], base.values_slice());
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
