#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::{evaluate_mixed, newton_tolerance};
  use approx::assert_abs_diff_eq;
  use ndarray::array;
  use num::clamp;
  use rstest::rstest;

  use super::super::test_newton_convergence_support::tests::make_dense_contribution;

  #[test]
  fn test_newton_tolerance_positive_branch_uses_relative() {
    // For branch_length >> NEWTON_ABS_TOL, tolerance is proportional
    let tol = newton_tolerance(0.1);
    assert_abs_diff_eq!(tol, 0.001 * 0.1, epsilon = 1e-15);
  }

  #[test]
  fn test_newton_tolerance_zero_branch_uses_absolute_floor() {
    // At branch_length=0, the relative component is 0 and the absolute
    // floor prevents degeneration
    let tol = newton_tolerance(0.0);
    assert_abs_diff_eq!(tol, 1e-8, epsilon = 1e-15);
    assert!(tol > 0.0, "Tolerance must be positive at zero branch length");
  }

  #[test]
  fn test_newton_tolerance_small_branch_uses_absolute_floor() {
    // For very small branches where 0.001 * bl < 1e-8, the floor dominates
    let tol = newton_tolerance(1e-7);
    assert_abs_diff_eq!(tol, 1e-8, epsilon = 1e-15);
  }

  #[test]
  fn test_newton_tolerance_crossover() {
    // At branch_length = 1e-5, relative = 1e-8 = absolute floor (crossover)
    let tol = newton_tolerance(1e-5);
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
  fn test_newton_tolerance_always_positive(#[case] bl: f64) {
    let tol = newton_tolerance(bl);
    assert!(tol > 0.0, "Tolerance must be positive for branch_length={bl}");
    assert!(tol.is_finite(), "Tolerance must be finite for branch_length={bl}");
  }

  #[test]
  fn test_newton_zero_start_converges_efficiently() {
    // Starting from branch_length=0, the absolute tolerance floor (NEWTON_ABS_TOL)
    // prevents degeneration to zero tolerance and allows early convergence.
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let mut branch_length = 0.0;
    let metrics = evaluate_mixed(&contributions, branch_length);
    let max_iter = 10;
    let mut n_iter = 0;

    if metrics.second_derivative < 0.0 {
      let mut new_branch_length =
        (branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length)).max(0.0);

      while (new_branch_length - branch_length).abs() > newton_tolerance(branch_length) && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length);
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

    // Must converge before exhausting all iterations
    assert!(
      n_iter < max_iter,
      "Newton loop from zero exhausted all {max_iter} iterations (n_iter={n_iter})"
    );
    assert!(branch_length >= 0.0);
    assert!(branch_length.is_finite());
  }
}
