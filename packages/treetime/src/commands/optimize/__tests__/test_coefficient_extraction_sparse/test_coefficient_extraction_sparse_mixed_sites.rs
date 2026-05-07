#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
  use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_mixed_variable_and_fixed_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    // Variable sites (multiplicity 1)
    let variable_coeffs = array![0.3, 0.3, 0.2, 0.2];
    let variable_site = SiteContribution {
      multiplicity: 1.0,
      coefficients: variable_coeffs.clone(),
    };

    // Fixed sites (multiplicity > 1)
    let fixed_coeffs = array![0.7, 0.1, 0.1, 0.1];
    let fixed_site = SiteContribution {
      multiplicity: 100.0,
      coefficients: fixed_coeffs.clone(),
    };

    let contribution = PartitionContribution {
      site_contributions: vec![variable_site, fixed_site],
      gtr: gtr.clone(),
    };

    let branch_length = 0.1;
    let metrics = evaluate_sparse_contribution(&contribution, branch_length);

    // Compute expected log-LH as sum of individual site contributions
    let variable_contribution = PartitionContribution {
      site_contributions: vec![SiteContribution {
        multiplicity: 1.0,
        coefficients: variable_coeffs,
      }],
      gtr: gtr.clone(),
    };
    let fixed_contribution = PartitionContribution {
      site_contributions: vec![SiteContribution {
        multiplicity: 100.0,
        coefficients: fixed_coeffs,
      }],
      gtr,
    };

    let variable_metrics = evaluate_sparse_contribution(&variable_contribution, branch_length);
    let fixed_metrics = evaluate_sparse_contribution(&fixed_contribution, branch_length);

    // Total log-LH equals sum of individual contributions
    let expected_log_lh = variable_metrics.log_lh + fixed_metrics.log_lh;
    pretty_assert_ulps_eq!(metrics.log_lh, expected_log_lh, max_ulps = 100);

    // Total derivative equals sum of individual derivatives
    let expected_derivative = variable_metrics.derivative + fixed_metrics.derivative;
    pretty_assert_ulps_eq!(metrics.derivative, expected_derivative, max_ulps = 100);
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
      gtr,
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
    pretty_assert_ulps_eq!(metrics.log_lh, expected, max_ulps = 100);
  }
}
