#[cfg(test)]
mod tests {
  use crate::optimize::likelihood::evaluate_mixed;
  use crate::optimize::method_newton::NEWTON_REL_TOL;
  use crate::optimize::method_newton::{newton_tolerance_log, newton_tolerance_sqrt, newton_tolerance_t};
  use approx::assert_abs_diff_eq;
  use ndarray::array;
  use num::clamp;
  use rstest::rstest;

  use super::super::test_newton_convergence_support::tests::make_dense_contribution;

  #[test]
  fn test_newton_tolerance_t_positive_branch_uses_relative() {
    let tol = newton_tolerance_t(0.1);
    assert_abs_diff_eq!(tol, 0.001 * 0.1, epsilon = 1e-15);
  }

  #[test]
  fn test_newton_tolerance_t_zero_branch_uses_absolute_floor() {
    let tol = newton_tolerance_t(0.0);
    assert_abs_diff_eq!(tol, 1e-8, epsilon = 1e-15);
    assert!(tol > 0.0, "Tolerance must be positive at zero branch length");
  }

  #[test]
  fn test_newton_tolerance_t_small_branch_uses_absolute_floor() {
    let tol = newton_tolerance_t(1e-7);
    assert_abs_diff_eq!(tol, 1e-8, epsilon = 1e-15);
  }

  #[test]
  fn test_newton_tolerance_t_crossover() {
    let tol = newton_tolerance_t(1e-5);
    assert_abs_diff_eq!(tol, 1e-8, epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(     0.0  )]
  #[case::tiny(     1e-15)]
  #[case::small(    1e-10)]
  #[case::at_floor( 1e-8 )]
  #[case::crossover(1e-5 )]
  #[case::moderate( 0.001)]
  #[case::typical(  0.1  )]
  #[case::one(      1.0  )]
  #[case::large(    10.0 )]
  #[trace]
  fn test_newton_tolerance_t_always_positive(#[case] bl: f64) {
    let tol = newton_tolerance_t(bl);
    assert!(tol > 0.0, "Tolerance must be positive for branch_length={bl}");
    assert!(tol.is_finite(), "Tolerance must be finite for branch_length={bl}");
  }

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

  #[test]
  fn test_newton_tolerance_log_is_constant() {
    let tol = newton_tolerance_log();
    assert_abs_diff_eq!(tol, NEWTON_REL_TOL.ln_1p(), epsilon = 1e-15);
    assert_abs_diff_eq!(tol, 0.001_f64.ln_1p(), epsilon = 1e-15);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small( 0.001)]
  #[case::medium(0.1  )]
  #[case::large( 1.0  )]
  #[trace]
  fn test_newton_tolerance_implied_relative_tolerance_consistent(#[case] t: f64) {
    let eta = NEWTON_REL_TOL;

    let implied_t = newton_tolerance_t(t) / t;

    let s = t.sqrt();
    let implied_sqrt = 2.0 * newton_tolerance_sqrt(s) / s;

    let implied_log = newton_tolerance_log();

    assert_abs_diff_eq!(implied_t, eta, epsilon = 1e-10);
    assert_abs_diff_eq!(implied_sqrt, eta, epsilon = 1e-10);
    assert_abs_diff_eq!(implied_log, eta, epsilon = 1e-6);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::moderate(1e-3 )]
  #[case::small(   1e-6 )]
  #[case::tiny(    1e-12)]
  #[trace]
  fn test_newton_tolerance_log_does_not_grow_with_short_branches(#[case] t: f64) {
    let _u = t.ln();
    let tol = newton_tolerance_log();
    assert!(tol < 0.002, "Log tolerance should be ≈0.001, got {tol} at t={t}");
    assert!(tol > 0.0009, "Log tolerance should be ≈0.001, got {tol} at t={t}");
  }

  #[test]
  fn test_newton_zero_start_converges_efficiently() {
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let mut branch_length = 0.0;
    let metrics = evaluate_mixed(&contributions, branch_length).expect("valid branch length");
    let max_iter = 10;
    let mut n_iter = 0;

    if metrics.second_derivative < 0.0 {
      let mut new_branch_length =
        (branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length)).max(0.0);

      while (new_branch_length - branch_length).abs() > newton_tolerance_t(branch_length) && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length).expect("valid branch length");
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = (branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            ))
          .max(0.0);
        } else {
          break;
        }
        n_iter += 1;
      }
      branch_length = new_branch_length;
    }

    assert!(
      n_iter < max_iter,
      "Newton loop from zero exhausted all {max_iter} iterations (n_iter={n_iter})"
    );
    assert!(branch_length >= 0.0);
    assert!(branch_length.is_finite());
  }
}
