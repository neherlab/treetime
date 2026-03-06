#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
  use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_evaluate_sparse_derivative_matches_numerical() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 5.0,
      coefficients: array![0.4, 0.3, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    // Verify analytical first derivative against numerical approximation
    let h = 1e-6;
    for &branch_length in &[0.01, 0.1, 0.5, 1.0] {
      let metrics = evaluate_sparse_contribution(&contribution, branch_length);

      // Numerical first derivative: (f(t+h) - f(t-h)) / 2h
      let metrics_plus = evaluate_sparse_contribution(&contribution, branch_length + h);
      let metrics_minus = evaluate_sparse_contribution(&contribution, branch_length - h);
      let numerical_d1 = (metrics_plus.log_lh - metrics_minus.log_lh) / (2.0 * h);

      assert_ulps_eq!(metrics.derivative, numerical_d1, epsilon = 1e-4);

      // Second derivative should be finite and negative (log-likelihood is concave)
      assert!(
        metrics.second_derivative.is_finite(),
        "second_derivative should be finite at branch_length={branch_length}"
      );
      assert!(
        metrics.second_derivative < 0.0,
        "second_derivative should be negative (concave) at branch_length={branch_length}"
      );
    }
  }

  #[test]
  fn test_evaluate_sparse_derivative_scales_with_multiplicity() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let coefficients = array![0.4, 0.3, 0.2, 0.1];

    let site1 = SiteContribution {
      multiplicity: 1.0,
      coefficients: coefficients.clone(),
    };
    let contribution1 = PartitionContribution {
      site_contributions: vec![site1],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let site3 = SiteContribution {
      multiplicity: 3.0,
      coefficients,
    };
    let contribution3 = PartitionContribution {
      site_contributions: vec![site3],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics1 = evaluate_sparse_contribution(&contribution1, 0.1);
    let metrics3 = evaluate_sparse_contribution(&contribution3, 0.1);

    // Derivative should scale with multiplicity
    assert_ulps_eq!(metrics3.derivative, 3.0 * metrics1.derivative, max_ulps = 100);
  }
}
