#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::optimize_dense::get_coefficients;
  use crate::pretty_assert_ulps_eq;
  use ndarray::array;

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  #[test]
  fn test_coefficients_use_eigenvector_decomposition() {
    // Verify that coefficients are computed via eigenvector decomposition:
    // k_c = (child · v)_c * (parent · v_inv^T)_c
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.4, 0.3, 0.2, 0.1]];
    let child = array![[0.1, 0.2, 0.3, 0.4]];
    let msg_to_parent = make_dense_seq_dis(parent.clone());
    let msg_to_child = make_dense_seq_dis(child.clone());

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Manually compute expected coefficients
    let child_v = child.dot(&gtr.v);
    let parent_v_inv_t = parent.dot(&gtr.v_inv.t());
    let expected = child_v * parent_v_inv_t;

    for i in 0..4 {
      pretty_assert_ulps_eq!(contribution.coefficients[[0, i]], expected[[0, i]], max_ulps = 10);
    }
  }

  #[test]
  fn test_coefficients_preserve_gtr_reference() {
    // PartitionContribution should store the GTR for later evaluation
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.25, 0.25, 0.25, 0.25]];
    let child = array![[0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    // Stored GTR should have same eigenvalues
    for i in 0..4 {
      pretty_assert_ulps_eq!(contribution.gtr.eigvals[i], gtr.eigvals[i], max_ulps = 10);
    }
  }
}
