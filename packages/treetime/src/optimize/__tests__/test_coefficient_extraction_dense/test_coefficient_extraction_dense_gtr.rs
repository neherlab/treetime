#[cfg(test)]
mod tests {
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::partition::optimize_dense::get_coefficients;
  use ndarray::array;
  use treetime_utils::pretty_assert_ulps_eq;

  use super::super::test_coefficient_extraction_dense_support::tests::make_dense_seq_dis;

  /// Coefficients are computed via eigenvector decomposition:
  /// k_c = (child . v)_c * (parent . v_inv^T)_c
  #[test]
  fn test_coefficients_use_eigenvector_decomposition() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.4, 0.3, 0.2, 0.1]];
    let child = array![[0.1, 0.2, 0.3, 0.4]];
    let msg_to_parent = make_dense_seq_dis(parent.clone());
    let msg_to_child = make_dense_seq_dis(child.clone());

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    let child_v = child.dot(&gtr.v);
    let parent_v_inv_t = parent.dot(&gtr.v_inv.t());
    let expected = child_v * parent_v_inv_t;

    pretty_assert_ulps_eq!(contribution.coefficients, expected, max_ulps = 10);
  }

  /// PartitionContribution should store the GTR for later evaluation.
  #[test]
  fn test_coefficients_preserve_gtr_reference() {
    let gtr = jc69(JC69Params::default()).expect("JC69 creation failed");

    let parent = array![[0.25, 0.25, 0.25, 0.25]];
    let child = array![[0.25, 0.25, 0.25, 0.25]];
    let msg_to_parent = make_dense_seq_dis(parent);
    let msg_to_child = make_dense_seq_dis(child);

    let contribution = get_coefficients(&msg_to_parent, &msg_to_child, &gtr);

    pretty_assert_ulps_eq!(contribution.gtr.eigvals, gtr.eigvals, max_ulps = 10);
  }
}
