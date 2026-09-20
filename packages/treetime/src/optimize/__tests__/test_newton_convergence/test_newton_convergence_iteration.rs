#[cfg(test)]
mod tests {
  use crate::optimize::likelihood::evaluate_with_indels;
  use crate::optimize::method_newton::newton_tolerance_t;
  use ndarray::array;
  use num::clamp;

  use super::super::test_newton_convergence_support::tests::make_dense_contribution;

  #[test]
  fn test_newton_iteration_converges_within_bounds() {
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let indel_count = 2_usize;
    let indel_rate = 10.0;
    let min_bl = 1e-5;
    let initial_bl = 0.1;
    let mut branch_length = initial_bl;
    let metrics =
      evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length).expect("valid branch length");
    let max_iter = 10;
    let mut n_iter = 0;

    assert!(
      metrics.second_derivative < 0.0,
      "precondition: second_derivative at branch_length=0.1 must be negative for Newton to run, got {}",
      metrics.second_derivative
    );
    let initial_lh = metrics.log_lh.value();

    let mut new_branch_length =
      (branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length)).max(min_bl);

    while (new_branch_length - branch_length).abs() > newton_tolerance_t(branch_length) && n_iter < max_iter {
      let new_metrics =
        evaluate_with_indels(&contributions, indel_count, indel_rate, new_branch_length).expect("valid branch length");
      if new_metrics.second_derivative < 0.0 {
        branch_length = new_branch_length;
        new_branch_length = (branch_length
          - clamp(
            new_metrics.derivative / new_metrics.second_derivative,
            -1.0,
            branch_length,
          ))
        .max(min_bl);
      } else {
        break;
      }
      n_iter += 1;
    }
    branch_length = new_branch_length;

    assert!(n_iter > 0, "expected at least one Newton iteration to run");
    assert!(n_iter <= max_iter);
    assert!(branch_length >= 0.0);
    let final_metrics =
      evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length).expect("valid branch length");
    assert!(
      final_metrics.log_lh.value() >= initial_lh - 1e-10,
      "Newton degraded log_lh from {initial_lh} to {final}",
      final = final_metrics.log_lh.value(),
    );
  }

  #[test]
  fn test_newton_iteration_respects_max_iter() {
    let coefficients = array![[0.5, 0.3, 0.1, 0.1], [0.4, 0.4, 0.1, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let indel_count = 2_usize;
    let indel_rate = 10.0;
    let min_bl = 1e-5;
    let mut branch_length = 0.5;
    let max_iter = 3;
    let mut n_iter = 0;

    let metrics =
      evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length).expect("valid branch length");
    assert!(
      metrics.second_derivative < 0.0,
      "precondition: chosen inputs must yield negative curvature, got {}",
      metrics.second_derivative
    );

    let mut new_branch_length =
      (branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length)).max(min_bl);

    let tight_tol = 1e-15;
    while (new_branch_length - branch_length).abs() > tight_tol && n_iter < max_iter {
      let new_metrics =
        evaluate_with_indels(&contributions, indel_count, indel_rate, new_branch_length).expect("valid branch length");
      if new_metrics.second_derivative < 0.0 {
        branch_length = new_branch_length;
        new_branch_length = (branch_length
          - clamp(
            new_metrics.derivative / new_metrics.second_derivative,
            -1.0,
            branch_length,
          ))
        .max(min_bl);
      } else {
        break;
      }
      n_iter += 1;
    }

    assert_eq!(
      max_iter, n_iter,
      "expected the iteration cap to fire, but loop exited early"
    );
  }
}
