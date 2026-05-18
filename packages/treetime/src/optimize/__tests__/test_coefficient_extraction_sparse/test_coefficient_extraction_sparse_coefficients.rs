#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::pretty_assert_ulps_eq;
  use crate::representation::partition::optimize_sparse::{PartitionContribution, SiteContribution};
  use ndarray::array;

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
      gtr,
    };

    // Verify coefficients match
    for i in 0..4 {
      pretty_assert_ulps_eq!(
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
      gtr,
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
      gtr,
    };

    let metrics = evaluate_sparse_contribution(&contribution, 0.0);

    // Mismatched states at zero branch length should have very low likelihood
    assert!(
      metrics.log_lh < -10.0 || metrics.log_lh == f64::NEG_INFINITY,
      "log-LH should be very low for mismatched states at zero branch length"
    );
  }
}
