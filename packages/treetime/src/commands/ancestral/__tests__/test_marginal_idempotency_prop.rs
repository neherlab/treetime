#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::arb_marginal_input_small;
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::representation::payload::ancestral::GraphAncestral;
  use proptest::prelude::*;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::prop_assert_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    /// Property test: `update_marginal` is a fixed-point operation on dense partitions
    /// across randomly generated trees, alignments, and GTR parameters.
    ///
    /// Felsenstein's pruning algorithm computes exact marginal likelihoods. Re-running
    /// it on the same data must yield the same log-likelihood. Random inputs exercise
    /// diverse tree topologies, branch lengths, sequence compositions, equilibrium
    /// frequencies, and substitution rate matrices to ensure the invariant holds
    /// regardless of parameterization.
    ///
    /// Invariant: `L(update_marginal(P)) == L(update_marginal(update_marginal(P)))`.
    ///
    /// Companion example test: `test_marginal_idempotency_example_dense`.
    #[test]
    fn test_prop_marginal_idempotency_dense(input in arb_marginal_input_small()) {
      let graph: GraphAncestral = nwk_read_str(&input.newick).unwrap();
      let (_, partitions) = run_dense_marginal(&input).unwrap();

      let log_lh_first = update_marginal(&graph, &partitions).unwrap();
      let log_lh_second = update_marginal(&graph, &partitions).unwrap();

      prop_assert_abs_diff_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);
    }

    /// Property test: `update_marginal` is a fixed-point operation on sparse partitions
    /// across randomly generated trees, alignments, and GTR parameters.
    ///
    /// The sparse representation compresses sequences by storing only variable positions
    /// and grouping fixed characters. This test verifies that the compression-aware
    /// marginal passes are idempotent, exercising different accumulation logic than the
    /// dense variant.
    ///
    /// Invariant: `L(update_marginal(P)) == L(update_marginal(update_marginal(P)))`.
    ///
    /// Companion example test: `test_marginal_idempotency_example_sparse`.
    #[test]
    fn test_prop_marginal_idempotency_sparse(input in arb_marginal_input_small()) {
      let graph: GraphAncestral = nwk_read_str(&input.newick).unwrap();
      let (_, partitions) = run_sparse_marginal(&input).unwrap();

      let log_lh_first = update_marginal(&graph, &partitions).unwrap();
      let log_lh_second = update_marginal(&graph, &partitions).unwrap();

      prop_assert_abs_diff_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);
    }
  }
}
