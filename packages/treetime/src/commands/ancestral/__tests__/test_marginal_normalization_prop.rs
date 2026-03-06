#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::{arb_marginal_input, arb_marginal_input_small};
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::representation::payload::sparse::MarginalSparseSeqDistribution;
  use ndarray::Array2;
  use proptest::prelude::*;

  /// Assert that every row of a dense profile matrix is a valid probability distribution.
  ///
  /// Checks three conditions per row:
  /// - Row sum equals 1.0 within epsilon = 1e-8 (abs-diff).
  /// - All values are finite (no NaN or Inf).
  /// - All values >= -1e-14 (allows roundoff-level negative values from
  ///   element-wise division in the outgroup message computation).
  ///
  /// Returns a proptest-compatible `TestCaseError` on failure to support shrinking.
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

  /// Assert that a sparse marginal profile is a valid probability distribution.
  ///
  /// The sparse representation separates positions into two groups:
  /// - Variable positions: sites that differ across the subtree, stored individually.
  /// - Fixed-character distributions: conditional posteriors for positions sharing
  ///   the same consensus state from Fitch compression, stored per character key.
  ///
  /// Checks the following for each component:
  /// - `log_lh` is finite.
  /// - Every variable-position distribution sums to 1.0 within epsilon = 1e-8 (abs-diff).
  /// - Every fixed-character distribution sums to 1.0 within epsilon = 1e-8 (abs-diff).
  /// - All values are finite and >= -1e-14 (allows roundoff-level negative values).
  ///
  /// Returns a proptest-compatible `TestCaseError` on failure to support shrinking.
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

    /// Property test: every node posterior and edge outgroup message in a dense marginal
    /// reconstruction is a valid probability distribution, across random trees and GTR models.
    ///
    /// Felsenstein's two-pass pruning algorithm (sum-product belief propagation on trees)
    /// guarantees that node posteriors and edge messages are normalizable at every step:
    /// - Backward pass (leaves to root): partial likelihoods are products of non-negative
    ///   transition-weighted child messages, normalized per row.
    /// - Forward pass (root to leaves): node posteriors combine ingroup and outgroup
    ///   information, normalized per row. Edge outgroup messages (`msg_to_child`) are
    ///   computed as node_posterior / msg_from_child, then normalized.
    ///
    /// Verifies that the total log-likelihood L = sum_s ln(P(D_s)) is finite and
    /// non-positive (each site likelihood P(D_s) is in (0, 1], so ln(P(D_s)) <= 0).
    ///
    /// Uses small random inputs (3-4 taxa, 3-10 positions) for thorough shrinking.
    ///
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

    /// Property test: every node posterior and edge outgroup message in a sparse marginal
    /// reconstruction is a valid probability distribution, across random trees and GTR models.
    ///
    /// The sparse representation splits positions into variable (stored individually) and
    /// fixed (grouped by consensus character from Fitch compression). The forward pass
    /// explicitly normalizes both variable-position and fixed-character distributions,
    /// so both must sum to 1. Verifies that the total log-likelihood is finite and
    /// non-positive (same invariant as the dense variant).
    ///
    /// Uses small random inputs (3-4 taxa, 3-10 positions) for thorough shrinking.
    ///
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

    /// Property test: dense marginal log-likelihood is finite and non-positive across
    /// larger random inputs (3-6 taxa, 5-20 positions).
    ///
    /// The total log-likelihood L = sum_s ln(P(D_s)) is non-positive because each
    /// site likelihood P(D_s) is a probability in (0, 1]. Larger inputs than the
    /// full-profile normalization tests stress the logsumexp accumulation path for
    /// numerical stability (underflow from multiplying many small partial likelihoods,
    /// cancellation in log-space sums) without the cost of inspecting every profile row.
    ///
    /// Companion example test: `test_marginal_normalization_example_dense`.
    #[test]
    fn test_prop_marginal_normalization_dense_log_lh_finite(input in arb_marginal_input()) {
      let (log_lh, _) = run_dense_marginal(&input).unwrap();
      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    }

    /// Property test: sparse marginal log-likelihood is finite and non-positive across
    /// larger random inputs (3-6 taxa, 5-20 positions).
    ///
    /// Same log-likelihood invariant as the dense variant (L <= 0), exercised on the
    /// sparse representation. The sparse path accumulates log-likelihood separately
    /// over variable positions (individual per-site contributions) and fixed positions
    /// (per-character contributions weighted by Fitch compression counts). Larger
    /// inputs increase both the number of variable sites and the compression multipliers,
    /// stressing this split accumulation for numerical stability.
    ///
    /// Companion example test: `test_marginal_normalization_example_sparse`.
    #[test]
    fn test_prop_marginal_normalization_sparse_log_lh_finite(input in arb_marginal_input()) {
      let (log_lh, _) = run_sparse_marginal(&input).unwrap();
      prop_assert!(log_lh.is_finite(), "Log-likelihood non-finite: {log_lh}");
      prop_assert!(log_lh <= 0.0, "Log-likelihood should be non-positive: {log_lh}");
    }
  }
}
