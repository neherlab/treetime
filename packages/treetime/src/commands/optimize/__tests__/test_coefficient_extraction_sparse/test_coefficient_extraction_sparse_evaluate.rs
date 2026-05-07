#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_sparse::{PartitionContribution, SiteContribution};
  use crate::commands::optimize::optimize_sparse_eval::evaluate_sparse_contribution;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

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
      gtr,
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // At branch_length=0, exp(λt) = 1 for all eigenvalues
    // log-LH = ln(sum of coefficients) = ln(1.0) = 0.0
    let coeff_sum: f64 = 0.5 + 0.2 + 0.2 + 0.1;
    pretty_assert_ulps_eq!(metrics.log_lh, coeff_sum.ln(), max_ulps = 100);
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
      gtr,
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // log-LH = multiplicity * ln(sum of coefficients) = 10 * ln(1.0) = 0.0
    let coeff_sum: f64 = 0.5 + 0.2 + 0.2 + 0.1;
    let expected_log_lh = 10.0 * coeff_sum.ln();
    pretty_assert_ulps_eq!(metrics.log_lh, expected_log_lh, max_ulps = 100);
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
      gtr: gtr.clone(),
    };

    // Multiplicity 5
    let site5 = SiteContribution {
      multiplicity: 5.0,
      coefficients,
    };
    let contribution5 = PartitionContribution {
      site_contributions: vec![site5],
      gtr,
    };

    let metrics1 = evaluate_sparse_contribution(&contribution1, 0.1);
    let metrics5 = evaluate_sparse_contribution(&contribution5, 0.1);

    // log-LH should scale with multiplicity
    pretty_assert_ulps_eq!(metrics5.log_lh, 5.0 * metrics1.log_lh, max_ulps = 100);
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
      gtr: gtr.clone(),
    };
    let contribution_b = PartitionContribution {
      site_contributions: vec![SiteContribution {
        multiplicity: 1.0,
        coefficients: coefficients_b.clone(),
      }],
      gtr: gtr.clone(),
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
      gtr,
    };
    let metrics_both = evaluate_sparse_contribution(&contribution_both, 0.1);

    // log-LH should be sum of individual contributions
    pretty_assert_ulps_eq!(metrics_both.log_lh, metrics_a.log_lh + metrics_b.log_lh, max_ulps = 100);
  }
}
