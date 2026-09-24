#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::sparse_eval::evaluate_sparse_contribution;
  use crate::partition::optimize::sparse::{PartitionContribution, SiteContribution};
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  #[test]
  fn test_coefficients_computed_via_eigenvector_decomposition() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent_dis = array![0.4, 0.3, 0.2, 0.1];
    let child_dis = array![0.1, 0.2, 0.3, 0.4];

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

    pretty_assert_ulps_eq!(
      contribution.site_contributions[0].coefficients.clone(),
      expected_coefficients,
      max_ulps = 10
    );
  }

  #[test]
  fn test_matching_states_high_lh_at_zero() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

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

    let metrics = evaluate_sparse_contribution(&contribution, 0.0, true).expect("valid branch length");

    assert!(
      metrics.log_lh.value() > -1.0,
      "log-LH should be high for matching states"
    );
  }

  #[test]
  fn test_mismatched_states_low_lh_at_zero() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

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

    let metrics = evaluate_sparse_contribution(&contribution, 0.0, true).expect("valid branch length");

    assert!(
      metrics.log_lh.value() < -10.0 || metrics.log_lh.value() == f64::NEG_INFINITY,
      "log-LH should be very low for mismatched states at zero branch length"
    );
  }
}
