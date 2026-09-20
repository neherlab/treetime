#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::{arb_marginal_input, arb_marginal_input_small};
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::partition::storage::sparse::SparseSeqDistribution;
  use ndarray::Array2;
  use proptest::prelude::*;
  use treetime_utils::{prop_assert_array_finite, prop_assert_array_nonneg};

  fn assert_dense_rows_normalized(dis: &Array2<f64>) -> Result<(), TestCaseError> {
    prop_assert_array_finite!(dis);
    prop_assert_array_nonneg!(dis, epsilon = 1e-14);
    for (row_idx, row) in dis.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      if !approx::abs_diff_eq!(sum, 1.0, epsilon = 1e-8) {
        return Err(TestCaseError::fail(format!(
          "row {row_idx}: sum = {sum}, expected 1.0 (epsilon = 1e-8)"
        )));
      }
    }
    Ok(())
  }

  fn assert_sparse_profile_normalized(profile: &SparseSeqDistribution) -> Result<(), TestCaseError> {
    prop_assert!(
      profile.log_lh.value().is_finite(),
      "Profile log_lh non-finite: {}",
      profile.log_lh.value()
    );

    for (pos, var_pos) in &profile.variable {
      prop_assert_array_finite!(var_pos.dis);
      prop_assert_array_nonneg!(var_pos.dis, epsilon = 1e-14);
      let sum: f64 = var_pos.dis.sum();
      if !approx::abs_diff_eq!(sum, 1.0, epsilon = 1e-8) {
        return Err(TestCaseError::fail(format!(
          "variable position {pos}: sum = {sum}, expected 1.0 (epsilon = 1e-8)"
        )));
      }
    }

    for (char_key, fixed_dis) in &profile.fixed {
      prop_assert_array_finite!(fixed_dis);
      prop_assert_array_nonneg!(fixed_dis, epsilon = 1e-14);
      let sum: f64 = fixed_dis.sum();
      if !approx::abs_diff_eq!(sum, 1.0, epsilon = 1e-8) {
        return Err(TestCaseError::fail(format!(
          "fixed distribution for {char_key:?}: sum = {sum}, expected 1.0 (epsilon = 1e-8)"
        )));
      }
    }
    Ok(())
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn test_prop_marginal_normalization_dense(input in arb_marginal_input_small()) {
      let (log_lh, partitions) = run_dense_marginal(&input).unwrap();

      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be <= 0: {log_lh}");

      let partition = &partitions;
      for node_data in partition.node_states.values() {
        if !node_data.profile.dis.is_empty() {
          assert_dense_rows_normalized(&node_data.profile.dis)?;
        }
      }
      for edge_data in partition.edges.forward.values() {
        if !edge_data.msg_to_child.dis.is_empty() {
          assert_dense_rows_normalized(&edge_data.msg_to_child.dis)?;
        }
      }
    }

    #[test]
    fn test_prop_marginal_normalization_sparse(input in arb_marginal_input_small()) {
      let (log_lh, partitions) = run_sparse_marginal(&input).unwrap();

      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be <= 0: {log_lh}");

      let partition = &partitions;
      for node_data in partition.node_states.values() {
        assert_sparse_profile_normalized(&node_data.profile)?;
      }
      for edge_data in partition.edges.forward.values() {
        assert_sparse_profile_normalized(&edge_data.msg_to_child)?;
      }
    }

    #[test]
    fn test_prop_marginal_normalization_dense_log_lh_finite(input in arb_marginal_input()) {
      let (log_lh, _) = run_dense_marginal(&input).unwrap();
      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    }

    #[test]
    fn test_prop_marginal_normalization_sparse_log_lh_finite(input in arb_marginal_input()) {
      let (log_lh, _) = run_sparse_marginal(&input).unwrap();
      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    }
  }
}
