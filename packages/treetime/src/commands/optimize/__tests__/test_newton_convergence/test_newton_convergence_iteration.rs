#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_unified::{evaluate_with_indels, newton_tolerance_t};
  use ndarray::array;
  use num::clamp;

  use super::super::test_newton_convergence_support::tests::make_dense_contribution;

  #[test]
  fn test_newton_iteration_converges_within_bounds() {
    // Simulate the Newton loop from optimize_unified.rs
    let coefficients = array![[0.9, 0.03, 0.03, 0.04], [0.03, 0.9, 0.03, 0.04], [0.1, 0.1, 0.7, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    // Indel contribution makes the second derivative negative even when the
    // substitution variance Var(lambda) is positive, so Newton actually runs.
    let indel_count = 2usize;
    let indel_rate = 10.0;
    // Floor matches min_branch_length_for_indels(2, one_mutation=0.001) so the
    // simulated loop preserves the Poisson indel domain (t > 0).
    let min_bl = 1e-5;
    let initial_bl = 0.1;
    let mut branch_length = initial_bl;
    let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
    let max_iter = 10;
    let mut n_iter = 0;

    // Precondition: the chosen inputs must have negative curvature at the
    // start, otherwise the test never enters the Newton loop and the rest
    // of the assertions become vacuous.
    assert!(
      metrics.second_derivative < 0.0,
      "precondition: second_derivative at branch_length=0.1 must be negative for Newton to run, got {}",
      metrics.second_derivative
    );
    let initial_lh = metrics.log_lh;

    let mut new_branch_length =
      (branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length)).max(min_bl);

    while (new_branch_length - branch_length).abs() > newton_tolerance_t(branch_length) && n_iter < max_iter {
      let new_metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, new_branch_length);
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

    // At least one Newton update ran (the loop body executed). Without this
    // the test would pass even if the loop guard rejected the very first
    // iteration.
    assert!(n_iter > 0, "expected at least one Newton iteration to run");
    // Standard envelope checks.
    assert!(n_iter <= max_iter);
    assert!(branch_length >= 0.0);
    // Newton must improve (or not degrade) the log-likelihood from the
    // starting point. This is the actual contract being tested.
    let final_metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
    assert!(
      final_metrics.log_lh >= initial_lh - 1e-10,
      "Newton degraded log_lh from {initial_lh} to {final}",
      final = final_metrics.log_lh,
    );
  }

  #[test]
  fn test_newton_iteration_respects_max_iter() {
    // Even with slow convergence, should stop at max_iter. Use a tiny
    // tolerance to force the loop to exhaust the iteration budget.
    let coefficients = array![[0.5, 0.3, 0.1, 0.1], [0.4, 0.4, 0.1, 0.1],];
    let contribution = make_dense_contribution(coefficients);
    let contributions = vec![contribution];

    let indel_count = 2usize;
    let indel_rate = 10.0;
    // Floor matches min_branch_length_for_indels(2, one_mutation=0.001) so the
    // simulated loop preserves the Poisson indel domain (t > 0).
    let min_bl = 1e-5;
    let mut branch_length = 0.5;
    let max_iter = 3;
    let mut n_iter = 0;

    let metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, branch_length);
    // Precondition: Newton must enter the loop, otherwise the max_iter
    // budget is never tested.
    assert!(
      metrics.second_derivative < 0.0,
      "precondition: chosen inputs must yield negative curvature, got {}",
      metrics.second_derivative
    );

    let mut new_branch_length =
      (branch_length - clamp(metrics.derivative / metrics.second_derivative, -1.0, branch_length)).max(min_bl);

    // Force the loop to run by using an impossibly tight per-step tolerance.
    let tight_tol = 1e-15;
    while (new_branch_length - branch_length).abs() > tight_tol && n_iter < max_iter {
      let new_metrics = evaluate_with_indels(&contributions, indel_count, indel_rate, new_branch_length);
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

    // The contract: the iteration count is bounded by max_iter even when
    // the per-step tolerance is impossibly tight.
    assert_eq!(
      max_iter, n_iter,
      "expected the iteration cap to fire, but loop exited early"
    );
  }
}
