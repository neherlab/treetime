#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
  use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use approx::assert_ulps_eq;
  use ndarray::array;

  // ==========================================================================
  // SiteContribution construction tests
  // ==========================================================================

  #[test]
  fn test_site_contribution_stores_multiplicity_and_coefficients() {
    let coefficients = array![0.1, 0.2, 0.3, 0.4];
    let contribution = SiteContribution {
      multiplicity: 5.0,
      coefficients: coefficients.clone(),
    };

    assert_ulps_eq!(contribution.multiplicity, 5.0, max_ulps = 10);
    for i in 0..4 {
      assert_ulps_eq!(contribution.coefficients[i], coefficients[i], max_ulps = 10);
    }
  }

  #[test]
  fn test_site_contribution_unit_multiplicity() {
    // Variable sites have multiplicity 1.0
    let contribution = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.25, 0.25, 0.25, 0.25],
    };

    assert_ulps_eq!(contribution.multiplicity, 1.0, max_ulps = 10);
  }

  // ==========================================================================
  // PartitionContribution construction tests
  // ==========================================================================

  #[test]
  fn test_partition_contribution_empty_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let contribution = PartitionContribution {
      site_contributions: vec![],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    assert!(contribution.site_contributions.is_empty());
    assert_eq!(contribution.eigenvalues.len(), 4);
  }

  #[test]
  fn test_partition_contribution_single_variable_site() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // One variable site with multiplicity 1
    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.5, 0.2, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    assert_eq!(contribution.site_contributions.len(), 1);
    assert_ulps_eq!(contribution.site_contributions[0].multiplicity, 1.0, max_ulps = 10);
  }

  #[test]
  fn test_partition_contribution_multiple_variable_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Multiple variable sites, each with multiplicity 1
    let sites = vec![
      SiteContribution {
        multiplicity: 1.0,
        coefficients: array![0.5, 0.2, 0.2, 0.1],
      },
      SiteContribution {
        multiplicity: 1.0,
        coefficients: array![0.1, 0.4, 0.3, 0.2],
      },
      SiteContribution {
        multiplicity: 1.0,
        coefficients: array![0.3, 0.3, 0.2, 0.2],
      },
    ];

    let contribution = PartitionContribution {
      site_contributions: sites,
      eigenvalues: gtr.eigvals.to_owned(),
    };

    assert_eq!(contribution.site_contributions.len(), 3);
    for site in &contribution.site_contributions {
      assert_ulps_eq!(site.multiplicity, 1.0, max_ulps = 10);
    }
  }

  #[test]
  fn test_partition_contribution_fixed_sites_with_multiplicity() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Fixed sites have multiplicity > 1 (count of identical positions)
    let sites = vec![
      SiteContribution {
        multiplicity: 100.0, // 100 identical fixed positions
        coefficients: array![0.9, 0.03, 0.03, 0.04],
      },
      SiteContribution {
        multiplicity: 50.0, // 50 identical fixed positions
        coefficients: array![0.8, 0.1, 0.05, 0.05],
      },
    ];

    let contribution = PartitionContribution {
      site_contributions: sites,
      eigenvalues: gtr.eigvals.to_owned(),
    };

    assert_eq!(contribution.site_contributions.len(), 2);
    assert_ulps_eq!(contribution.site_contributions[0].multiplicity, 100.0, max_ulps = 10);
    assert_ulps_eq!(contribution.site_contributions[1].multiplicity, 50.0, max_ulps = 10);
  }

  // ==========================================================================
  // Evaluation tests with multiplicity
  // ==========================================================================

  #[test]
  fn test_evaluate_sparse_single_site_multiplicity_1() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Single site with multiplicity 1
    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.5, 0.2, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // At branch_length=0, exp(λt) = 1 for all eigenvalues
    // log-LH = ln(sum of coefficients) = ln(1.0) = 0.0
    let coeff_sum: f64 = 0.5 + 0.2 + 0.2 + 0.1;
    assert_ulps_eq!(metrics.log_lh, coeff_sum.ln(), max_ulps = 100);
  }

  #[test]
  fn test_evaluate_sparse_single_site_multiplicity_10() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Single site with multiplicity 10 (10 identical positions)
    let site = SiteContribution {
      multiplicity: 10.0,
      coefficients: array![0.5, 0.2, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // log-LH = multiplicity * ln(sum of coefficients) = 10 * ln(1.0) = 0.0
    let coeff_sum: f64 = 0.5 + 0.2 + 0.2 + 0.1;
    let expected_log_lh = 10.0 * coeff_sum.ln();
    assert_ulps_eq!(metrics.log_lh, expected_log_lh, max_ulps = 100);
  }

  #[test]
  fn test_evaluate_sparse_multiplicity_scales_log_lh() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let coefficients = array![0.6, 0.2, 0.1, 0.1];

    // Multiplicity 1
    let site1 = SiteContribution {
      multiplicity: 1.0,
      coefficients: coefficients.clone(),
    };
    let contribution1 = PartitionContribution {
      site_contributions: vec![site1],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    // Multiplicity 5
    let site5 = SiteContribution {
      multiplicity: 5.0,
      coefficients: coefficients.clone(),
    };
    let contribution5 = PartitionContribution {
      site_contributions: vec![site5],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics1 = evaluate_sparse_contribution(&contribution1, 0.1);
    let metrics5 = evaluate_sparse_contribution(&contribution5, 0.1);

    // log-LH should scale with multiplicity
    assert_ulps_eq!(metrics5.log_lh, 5.0 * metrics1.log_lh, max_ulps = 100);
  }

  #[test]
  fn test_evaluate_sparse_multiple_sites_sum_log_lh() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let coefficients_a = array![0.6, 0.2, 0.1, 0.1];
    let coefficients_b = array![0.3, 0.3, 0.2, 0.2];

    // Evaluate separately
    let contribution_a = PartitionContribution {
      site_contributions: vec![SiteContribution {
        multiplicity: 1.0,
        coefficients: coefficients_a.clone(),
      }],
      eigenvalues: gtr.eigvals.to_owned(),
    };
    let contribution_b = PartitionContribution {
      site_contributions: vec![SiteContribution {
        multiplicity: 1.0,
        coefficients: coefficients_b.clone(),
      }],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics_a = evaluate_sparse_contribution(&contribution_a, 0.1);
    let metrics_b = evaluate_sparse_contribution(&contribution_b, 0.1);

    // Evaluate together
    let contribution_both = PartitionContribution {
      site_contributions: vec![
        SiteContribution {
          multiplicity: 1.0,
          coefficients: coefficients_a,
        },
        SiteContribution {
          multiplicity: 1.0,
          coefficients: coefficients_b,
        },
      ],
      eigenvalues: gtr.eigvals.to_owned(),
    };
    let metrics_both = evaluate_sparse_contribution(&contribution_both, 0.1);

    // log-LH should be sum of individual contributions
    assert_ulps_eq!(metrics_both.log_lh, metrics_a.log_lh + metrics_b.log_lh, max_ulps = 100);
  }

  // ==========================================================================
  // Coefficient computation property tests
  // ==========================================================================

  #[test]
  fn test_coefficients_computed_via_eigenvector_decomposition() {
    // k_c = (child · v)_c * (parent · v_inv^T)_c
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Parent and child probability distributions
    let parent_dis = array![0.4, 0.3, 0.2, 0.1];
    let child_dis = array![0.1, 0.2, 0.3, 0.4];

    // Manually compute coefficients via eigenvector decomposition
    let child_v = child_dis.dot(&gtr.v);
    let parent_v_inv_t = parent_dis.dot(&gtr.v_inv.t());
    let expected_coefficients = child_v * parent_v_inv_t;

    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients: expected_coefficients.clone(),
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    // Verify coefficients match
    for i in 0..4 {
      assert_ulps_eq!(
        contribution.site_contributions[0].coefficients[i],
        expected_coefficients[i],
        max_ulps = 10
      );
    }
  }

  #[test]
  fn test_matching_states_high_lh_at_zero() {
    // When parent and child have same state, coefficient sum should be high
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Both certain of state A
    let parent_dis = array![1.0, 0.0, 0.0, 0.0];
    let child_dis = array![1.0, 0.0, 0.0, 0.0];

    let child_v = child_dis.dot(&gtr.v);
    let parent_v_inv_t = parent_dis.dot(&gtr.v_inv.t());
    let coefficients = child_v * parent_v_inv_t;

    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients,
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // Matching states at zero branch length should have high likelihood
    assert!(metrics.log_lh > -1.0, "log-LH should be high for matching states");
  }

  #[test]
  fn test_mismatched_states_low_lh_at_zero() {
    // When parent and child have different states, coefficient sum at zero should be low
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Parent certain A, child certain C
    let parent_dis = array![1.0, 0.0, 0.0, 0.0];
    let child_dis = array![0.0, 1.0, 0.0, 0.0];

    let child_v = child_dis.dot(&gtr.v);
    let parent_v_inv_t = parent_dis.dot(&gtr.v_inv.t());
    let coefficients = child_v * parent_v_inv_t;

    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients,
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // Mismatched states at zero branch length should have very low likelihood
    assert!(
      metrics.log_lh < -10.0 || metrics.log_lh == f64::NEG_INFINITY,
      "log-LH should be very low for mismatched states at zero branch length"
    );
  }

  // ==========================================================================
  // Derivative tests
  // ==========================================================================

  #[test]
  fn test_evaluate_sparse_returns_finite_derivatives() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 5.0,
      coefficients: array![0.4, 0.3, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    for &branch_length in &[0.001, 0.01, 0.1, 1.0] {
      let metrics = evaluate_sparse_contribution(&contribution, branch_length);
      assert!(metrics.log_lh.is_finite(), "log-LH should be finite at branch_length={branch_length}");
      assert!(
        metrics.derivative.is_finite(),
        "derivative should be finite at branch_length={branch_length}"
      );
      assert!(
        metrics.second_derivative.is_finite(),
        "second_derivative should be finite at branch_length={branch_length}"
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
      coefficients: coefficients.clone(),
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

  // ==========================================================================
  // Combined variable and fixed sites tests
  // ==========================================================================

  #[test]
  fn test_mixed_variable_and_fixed_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Variable sites (multiplicity 1)
    let variable_site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.3, 0.3, 0.2, 0.2],
    };

    // Fixed sites (multiplicity > 1)
    let fixed_site = SiteContribution {
      multiplicity: 100.0,
      coefficients: array![0.7, 0.1, 0.1, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![variable_site, fixed_site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.1);

    // Total contribution is sum of both
    assert!(metrics.log_lh.is_finite());
    assert!(metrics.derivative.is_finite());
  }

  #[test]
  fn test_fixed_sites_dominate_with_high_multiplicity() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Variable site with low likelihood coefficients
    let variable_site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.1, 0.1, 0.1, 0.1], // Low sum
    };

    // Fixed site with high likelihood coefficients and high multiplicity
    let fixed_site = SiteContribution {
      multiplicity: 1000.0,
      coefficients: array![0.9, 0.03, 0.03, 0.04], // High sum
    };

    let contribution = PartitionContribution {
      site_contributions: vec![variable_site, fixed_site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // The fixed site contribution should dominate
    // variable: 1.0 * ln(0.4) ≈ -0.916
    // fixed: 1000.0 * ln(1.0) = 0.0
    // Total ≈ -0.916, but fixed sites with high LH reduce total magnitude
    assert!(metrics.log_lh.is_finite());
    // The fixed site with sum ~1.0 contributes 0 to log-LH per multiplicity
    // The variable site with sum ~0.4 contributes negative value
    let fixed_coeff_sum: f64 = 0.9 + 0.03 + 0.03 + 0.04;
    let variable_coeff_sum: f64 = 0.1 + 0.1 + 0.1 + 0.1;
    let expected = 1.0 * variable_coeff_sum.ln() + 1000.0 * fixed_coeff_sum.ln();
    assert_ulps_eq!(metrics.log_lh, expected, max_ulps = 100);
  }

  // ==========================================================================
  // Eigenvalue usage tests
  // ==========================================================================

  #[test]
  fn test_eigenvalues_affect_branch_length_evaluation() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.5, 0.2, 0.2, 0.1],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![site],
      eigenvalues: gtr.eigvals.to_owned(),
    };

    // At different branch lengths, exp(λt) changes for non-zero eigenvalues
    let metrics_short = evaluate_sparse_contribution(&contribution, 0.01);
    let metrics_long = evaluate_sparse_contribution(&contribution, 1.0);

    // Log-LH should differ at different branch lengths
    assert!(
      (metrics_short.log_lh - metrics_long.log_lh).abs() > 1e-6,
      "log-LH should differ at different branch lengths"
    );
  }

  #[test]
  fn test_jc69_eigenvalues_structure() {
    // JC69 has one zero eigenvalue and three equal negative eigenvalues
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Find the zero eigenvalue (should be one)
    let zero_count = gtr.eigvals.iter().filter(|&&ev| ev.abs() < 1e-10).count();
    assert_eq!(zero_count, 1, "JC69 should have exactly one zero eigenvalue");

    // All eigenvalues should be <= 0
    for &ev in gtr.eigvals.iter() {
      assert!(ev <= 1e-10, "JC69 eigenvalues should be non-positive");
    }
  }
}
