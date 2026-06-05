#[cfg(test)]
mod tests {
  use crate::piecewise_fn::PiecewiseFnBase;
  use ndarray::array;

  #[test]
  fn test_piecewise_fn_base_accessors() {
    let base = PiecewiseFnBase::new(array![1.0, 5.0, 10.0], array![0.0, 1.0, 2.0, 3.0]);
    assert_eq!(base.breakpoints(), &array![1.0, 5.0, 10.0]);
    assert_eq!(base.values(), &array![0.0, 1.0, 2.0, 3.0]);
    assert_eq!(base.breakpoints_slice(), &[1.0, 5.0, 10.0]);
    assert_eq!(base.values_slice(), &[0.0, 1.0, 2.0, 3.0]);
  }

  #[test]
  fn test_piecewise_fn_base_single_breakpoint() {
    let base = PiecewiseFnBase::new(array![5.0], array![0.0, 1.0]);
    assert_eq!(base.breakpoints().len(), 1);
    assert_eq!(base.values().len(), 2);
  }

  #[test]
  fn test_piecewise_fn_base_clone() {
    let base = PiecewiseFnBase::new(array![1.0, 2.0], array![10.0, 20.0]);
    let cloned = base.clone();
    assert_eq!(cloned.breakpoints(), base.breakpoints());
    assert_eq!(cloned.values(), base.values());
  }
}
