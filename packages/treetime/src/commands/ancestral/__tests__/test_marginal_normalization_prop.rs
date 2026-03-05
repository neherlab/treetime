#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::{arb_marginal_input, arb_marginal_input_small};
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::representation::payload::sparse::MarginalSparseSeqDistribution;
  use ndarray::Array2;
  use proptest::prelude::*;

  fn assert_dense_rows_normalized(dis: &Array2<f64>) -> Result<(), TestCaseError> {
    for (row_idx, row) in dis.rows().into_iter().enumerate() {
      let sum: f64 = row.sum();
      if !approx::abs_diff_eq!(sum, 1.0, epsilon = 1e-8) {
        return Err(TestCaseError::fail(format!(
          "row {row_idx}: sum = {sum}, expected 1.0 (epsilon = 1e-8)"
        )));
      }
      for (col_idx, &val) in row.iter().enumerate() {
        prop_assert!(val.is_finite(), "row {row_idx}, col {col_idx} non-finite: {val}");
        prop_assert!(val >= -1e-14, "row {row_idx}, col {col_idx} negative: {val}");
      }
    }
    Ok(())
  }

  fn assert_sparse_profile_normalized(profile: &MarginalSparseSeqDistribution) -> Result<(), TestCaseError> {
    prop_assert!(
      profile.log_lh.is_finite(),
      "Profile log_lh non-finite: {}",
      profile.log_lh
    );

    for (pos, var_pos) in &profile.variable {
      let sum: f64 = var_pos.dis.sum();
      if !approx::abs_diff_eq!(sum, 1.0, epsilon = 1e-8) {
        return Err(TestCaseError::fail(format!(
          "variable position {pos}: sum = {sum}, expected 1.0 (epsilon = 1e-8)"
        )));
      }
      for (idx, &val) in var_pos.dis.iter().enumerate() {
        prop_assert!(val.is_finite(), "variable pos {pos}, idx {idx} non-finite: {val}");
        prop_assert!(val >= -1e-14, "variable pos {pos}, idx {idx} negative: {val}");
      }
    }

    for (char_key, fixed_dis) in &profile.fixed {
      let sum: f64 = fixed_dis.sum();
      if !approx::abs_diff_eq!(sum, 1.0, epsilon = 1e-8) {
        return Err(TestCaseError::fail(format!(
          "fixed distribution for {char_key:?}: sum = {sum}, expected 1.0 (epsilon = 1e-8)"
        )));
      }
      for (idx, &val) in fixed_dis.iter().enumerate() {
        prop_assert!(val.is_finite(), "fixed {char_key:?}, idx {idx} non-finite: {val}");
        prop_assert!(val >= -1e-14, "fixed {char_key:?}, idx {idx} negative: {val}");
      }
    }
    Ok(())
  }

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    /// Companion example test: `test_marginal_normalization_example_dense`.
    #[test]
    fn test_prop_marginal_normalization_dense(input in arb_marginal_input_small()) {
      let (log_lh, partitions) = run_dense_marginal(&input).unwrap();

      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be <= 0: {log_lh}");

      let partition = partitions[0].read_arc();
      for node_data in partition.nodes.values() {
        if !node_data.profile.dis.is_empty() {
          assert_dense_rows_normalized(&node_data.profile.dis)?;
        }
      }
      for edge_data in partition.edges.values() {
        if !edge_data.msg_to_child.dis.is_empty() {
          assert_dense_rows_normalized(&edge_data.msg_to_child.dis)?;
        }
      }
    }

    /// Companion example test: `test_marginal_normalization_example_sparse`.
    #[test]
    fn test_prop_marginal_normalization_sparse(input in arb_marginal_input_small()) {
      let (log_lh, partitions) = run_sparse_marginal(&input).unwrap();

      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be <= 0: {log_lh}");

      let partition = partitions[0].read_arc();
      for node_data in partition.nodes.values() {
        assert_sparse_profile_normalized(&node_data.profile)?;
      }
      for edge_data in partition.edges.values() {
        assert_sparse_profile_normalized(&edge_data.msg_to_child)?;
      }
    }

    /// Companion example test: `test_marginal_normalization_example_dense`.
    #[test]
    fn test_prop_marginal_normalization_dense_log_lh_finite(input in arb_marginal_input()) {
      let (log_lh, _) = run_dense_marginal(&input).unwrap();
      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    }

    /// Companion example test: `test_marginal_normalization_example_sparse`.
    #[test]
    fn test_prop_marginal_normalization_sparse_log_lh_finite(input in arb_marginal_input()) {
      let (log_lh, _) = run_sparse_marginal(&input).unwrap();
      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    }
  }
}
