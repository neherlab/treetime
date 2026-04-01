#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_dense;
  use crate::commands::optimize::optimize_dense_eval::evaluate_dense_contribution;
  use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
  use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use ndarray::{Array2, array};
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::short( 0.01)]
  #[case::medium(0.1 )]
  #[case::long(  0.5 )]
  #[case::very_long(1.0)]
  #[trace]
  fn test_evaluate_sparse_first_derivative_matches_numerical(#[case] branch_length: f64) {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 5.0,
      coefficients: array![0.4, 0.3, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
      unimodal_branch_likelihood: gtr.unimodal_branch_likelihood,
    };

    let h = 1e-6;
    let metrics = evaluate_sparse_contribution(&contribution, branch_length);

    // Numerical first derivative: (f(t+h) - f(t-h)) / 2h
    let metrics_plus = evaluate_sparse_contribution(&contribution, branch_length + h);
    let metrics_minus = evaluate_sparse_contribution(&contribution, branch_length - h);
    let numerical_d1 = (metrics_plus.log_lh - metrics_minus.log_lh) / (2.0 * h);

    assert_ulps_eq!(metrics.derivative, numerical_d1, epsilon = 1e-4);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::short( 0.01)]
  #[case::medium(0.1 )]
  #[case::long(  0.5 )]
  #[case::very_long(1.0)]
  #[trace]
  fn test_evaluate_sparse_second_derivative_matches_numerical(#[case] branch_length: f64) {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 5.0,
      coefficients: array![0.4, 0.3, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
      unimodal_branch_likelihood: gtr.unimodal_branch_likelihood,
    };

    // Numerical second derivative via central difference of first derivative:
    // d2 ≈ (d1(t+h) - d1(t-h)) / 2h
    let h = 1e-5;
    let metrics = evaluate_sparse_contribution(&contribution, branch_length);
    let metrics_plus = evaluate_sparse_contribution(&contribution, branch_length + h);
    let metrics_minus = evaluate_sparse_contribution(&contribution, branch_length - h);
    let numerical_d2 = (metrics_plus.derivative - metrics_minus.derivative) / (2.0 * h);

    assert_ulps_eq!(metrics.second_derivative, numerical_d2, epsilon = 1e-4);
  }

  /// Second derivative with high multiplicity must match numerical approximation.
  /// Before the fix, multiplicity was squared in the Hessian, causing divergence
  /// proportional to m for large m.
  #[rustfmt::skip]
  #[rstest]
  #[case::short( 0.01)]
  #[case::medium(0.1 )]
  #[case::long(  0.5 )]
  #[case::very_long(1.0)]
  #[trace]
  fn test_evaluate_sparse_second_derivative_numerical_high_multiplicity(#[case] branch_length: f64) {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 500.0,
      coefficients: array![0.4, 0.3, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
      unimodal_branch_likelihood: gtr.unimodal_branch_likelihood,
    };

    let h = 1e-5;
    let metrics = evaluate_sparse_contribution(&contribution, branch_length);
    let metrics_plus = evaluate_sparse_contribution(&contribution, branch_length + h);
    let metrics_minus = evaluate_sparse_contribution(&contribution, branch_length - h);
    let numerical_d2 = (metrics_plus.derivative - metrics_minus.derivative) / (2.0 * h);

    assert_ulps_eq!(metrics.second_derivative, numerical_d2, epsilon = 1e-4);
  }

  #[test]
  fn test_evaluate_sparse_derivatives_scale_with_multiplicity() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let coefficients = array![0.4, 0.3, 0.2, 0.1];

    let site1 = SiteContribution {
      multiplicity: 1.0,
      coefficients: coefficients.clone(),
    };
    let contribution1 = PartitionContribution {
      site_contributions: vec![site1],
      eigenvalues: gtr.eigvals.to_owned(),
      unimodal_branch_likelihood: gtr.unimodal_branch_likelihood,
    };

    let site3 = SiteContribution {
      multiplicity: 3.0,
      coefficients,
    };
    let contribution3 = PartitionContribution {
      site_contributions: vec![site3],
      eigenvalues: gtr.eigvals.to_owned(),
      unimodal_branch_likelihood: gtr.unimodal_branch_likelihood,
    };

    let metrics1 = evaluate_sparse_contribution(&contribution1, 0.1);
    let metrics3 = evaluate_sparse_contribution(&contribution3, 0.1);

    // All three quantities scale linearly with multiplicity:
    // ℓ(t) = m * ln(L), ℓ'(t) = m * d1, ℓ''(t) = m * (d2 - d1^2)
    assert_ulps_eq!(metrics3.log_lh, 3.0 * metrics1.log_lh, max_ulps = 100);
    assert_ulps_eq!(metrics3.derivative, 3.0 * metrics1.derivative, max_ulps = 100);
    assert_ulps_eq!(
      metrics3.second_derivative,
      3.0 * metrics1.second_derivative,
      max_ulps = 100
    );
  }

  /// Sparse with multiplicity=1 must produce identical metrics to dense for the
  /// same coefficients and eigenvalues.
  #[rustfmt::skip]
  #[rstest]
  #[case::short( 0.01)]
  #[case::medium(0.1 )]
  #[case::long(  0.5 )]
  #[case::very_long(1.0)]
  #[trace]
  fn test_evaluate_sparse_matches_dense_multiplicity_1(#[case] branch_length: f64) {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let coefficients_a = array![0.4, 0.3, 0.2, 0.1];
    let coefficients_b = array![0.3, 0.3, 0.2, 0.2];

    // Sparse: two sites with multiplicity 1
    let sparse_contribution = PartitionContribution {
      site_contributions: vec![
        SiteContribution {
          multiplicity: 1.0,
          coefficients: coefficients_a.clone(),
        },
        SiteContribution {
          multiplicity: 1.0,
          coefficients: coefficients_b.clone(),
        },
      ],
      eigenvalues: gtr.eigvals.to_owned(),
      unimodal_branch_likelihood: gtr.unimodal_branch_likelihood,
    };

    // Dense: same two sites as rows of a 2D array
    let coefficients_2d = Array2::from_shape_vec((2, 4), {
      let mut v = coefficients_a.to_vec();
      v.extend(coefficients_b.to_vec());
      v
    })
    .expect("shape mismatch");

    let dense_contribution = optimize_dense::PartitionContribution::new(coefficients_2d, gtr);

    let sparse_metrics = evaluate_sparse_contribution(&sparse_contribution, branch_length);
    let dense_metrics = evaluate_dense_contribution(&dense_contribution, branch_length);

    assert_ulps_eq!(sparse_metrics.log_lh, dense_metrics.log_lh, max_ulps = 10);
    assert_ulps_eq!(sparse_metrics.derivative, dense_metrics.derivative, max_ulps = 10);
    assert_ulps_eq!(
      sparse_metrics.second_derivative,
      dense_metrics.second_derivative,
      max_ulps = 10
    );
  }
}
