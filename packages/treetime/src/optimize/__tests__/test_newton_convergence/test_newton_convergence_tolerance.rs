#[cfg(test)]
mod tests {
  use crate::optimize::method_newton::newton_tolerance_sqrt;
  use approx::assert_abs_diff_eq;

  #[test]
  fn test_newton_tolerance_sqrt_half_factor() {
    let s = 0.3;
    let tol = newton_tolerance_sqrt(s);
    assert_abs_diff_eq!(tol, 0.5 * 0.001 * 0.3, epsilon = 1e-15);
  }

  #[test]
  fn test_newton_tolerance_sqrt_zero_uses_floor() {
    let tol = newton_tolerance_sqrt(0.0);
    assert_abs_diff_eq!(tol, 1e-8, epsilon = 1e-15);
  }
}
