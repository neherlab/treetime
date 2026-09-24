#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::partition::optimize::sparse::{PartitionContribution, SiteContribution};
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

  #[test]
  fn test_mixed_variable_and_fixed_sites() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let variable_coeffs = array![0.3, 0.3, 0.2, 0.2];
    let variable_site = SiteContribution {
      multiplicity: 1.0,
      coefficients: variable_coeffs.clone(),
    };

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
    let metrics = evaluate_sparse_contribution(&contribution, branch_length, true).expect("valid branch length");

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

    let variable_metrics =
      evaluate_sparse_contribution(&variable_contribution, branch_length, true).expect("valid branch length");
    let fixed_metrics =
      evaluate_sparse_contribution(&fixed_contribution, branch_length, true).expect("valid branch length");

    let expected_log_lh = variable_metrics.log_lh.value() + fixed_metrics.log_lh.value();
    pretty_assert_ulps_eq!(metrics.log_lh.value(), expected_log_lh, max_ulps = 100);

    let expected_derivative = variable_metrics.derivative + fixed_metrics.derivative;
    pretty_assert_ulps_eq!(metrics.derivative, expected_derivative, max_ulps = 100);
  }

  #[test]
  fn test_fixed_sites_dominate_with_high_multiplicity() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let variable_site = SiteContribution {
      multiplicity: 1.0,
      coefficients: array![0.1, 0.1, 0.1, 0.1],
    };

    let fixed_site = SiteContribution {
      multiplicity: 1000.0,
      coefficients: array![0.9, 0.03, 0.03, 0.04],
    };

    let contribution = PartitionContribution {
      site_contributions: vec![variable_site, fixed_site],
      gtr,
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0, true).expect("valid branch length");

    assert!(metrics.log_lh.value().is_finite());
    let fixed_coeff_sum: f64 = 0.9 + 0.03 + 0.03 + 0.04;
    let variable_coeff_sum: f64 = 0.1 + 0.1 + 0.1 + 0.1;
    let expected = 1.0 * variable_coeff_sum.ln() + 1000.0 * fixed_coeff_sum.ln();
    pretty_assert_ulps_eq!(metrics.log_lh.value(), expected, max_ulps = 100);
  }
}
