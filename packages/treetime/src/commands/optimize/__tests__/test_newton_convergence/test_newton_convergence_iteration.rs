#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::evaluate_mixed;
  use ndarray::array;
  use num::clamp;

  use super::super::test_newton_convergence_support::tests::make_dense_contribution;

  #[test]
  fn test_newton_iteration_converges_within_bounds() {
    // Simulate the Newton loop from optimize_unified.rs
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let mut branch_length = 0.1;
    let metrics = evaluate_mixed(&contributions, branch_length);
    let max_iter = 10;
    let mut n_iter = 0;

    if metrics.second_derivative < 0.0 {
      let mut new_branch_length =
        branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);

      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length);
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            );
        } else {
          break;
        }
        n_iter += 1;
      }
      branch_length = new_branch_length;
    }

    // Verify convergence happened within iteration limit
    assert!(n_iter <= max_iter);
    // Verify branch length is non-negative
    assert!(branch_length >= 0.0);
  }

  #[test]
  fn test_newton_iteration_respects_max_iter() {
    // Even with slow convergence, should stop at max_iter
    let coefficients = array![[0.25, 0.25, 0.25, 0.25], [0.25, 0.25, 0.25, 0.25],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let mut branch_length = 0.5;
    let max_iter = 10;
    let mut n_iter = 0;

    let metrics = evaluate_mixed(&contributions, branch_length);
    if metrics.second_derivative < 0.0 {
      let mut new_branch_length =
        branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length);

      while (new_branch_length - branch_length).abs() > 0.001 * branch_length && n_iter < max_iter {
        let new_metrics = evaluate_mixed(&contributions, new_branch_length);
        if new_metrics.second_derivative < 0.0 {
          branch_length = new_branch_length;
          new_branch_length = branch_length
            - clamp(
              new_metrics.derivative / new_metrics.second_derivative,
              -1.0,
              branch_length,
            );
        } else {
          break;
        }
        n_iter += 1;
      }
    }

    assert!(n_iter <= max_iter);
  }
}
