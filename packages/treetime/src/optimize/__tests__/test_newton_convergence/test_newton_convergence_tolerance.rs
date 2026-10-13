#[cfg(test)]
mod tests {
  use crate::optimize::method_newton::{NEWTON_REL_TOL, newton_tolerance_sqrt};
  use approx::assert_abs_diff_eq;
  use rstest::rstest;

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

  #[rustfmt::skip]
  #[rstest]
  #[case::small( 0.001)]
  #[case::medium(0.1  )]
  #[case::large( 1.0  )]
  #[trace]
  fn test_newton_tolerance_sqrt_implies_relative_tolerance_in_t(#[case] t: f64) {
    let s = t.sqrt();
    let implied_relative_t = 2.0 * newton_tolerance_sqrt(s) / s;

    assert_abs_diff_eq!(NEWTON_REL_TOL, implied_relative_t, epsilon = 1e-10);
  }
}
