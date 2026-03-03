#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::arb_marginal_input_small;
  use crate::commands::ancestral::__tests__::prop_marginal_support::{run_dense_marginal, run_sparse_marginal};
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::representation::payload::ancestral::GraphAncestral;
  use proptest::prelude::*;
  use treetime_io::nwk::nwk_read_str;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    /// Companion example test: `test_marginal_idempotency_example_dense`.
    #[test]
    fn test_prop_marginal_idempotency_dense(input in arb_marginal_input_small()) {
      let graph: GraphAncestral = nwk_read_str(&input.newick).unwrap();
      let (_, partitions) = run_dense_marginal(&input).unwrap();

      let log_lh_first = update_marginal(&graph, &partitions).unwrap();
      let log_lh_second = update_marginal(&graph, &partitions).unwrap();

      prop_assert!(
        (log_lh_first - log_lh_second).abs() < 1e-10,
        "Idempotency violated: {log_lh_first} != {log_lh_second}"
      );
    }

    /// Companion example test: `test_marginal_idempotency_example_sparse`.
    #[test]
    fn test_prop_marginal_idempotency_sparse(input in arb_marginal_input_small()) {
      let graph: GraphAncestral = nwk_read_str(&input.newick).unwrap();
      let (_, partitions) = run_sparse_marginal(&input).unwrap();

      let log_lh_first = update_marginal(&graph, &partitions).unwrap();
      let log_lh_second = update_marginal(&graph, &partitions).unwrap();

      prop_assert!(
        (log_lh_first - log_lh_second).abs() < 1e-10,
        "Idempotency violated: {log_lh_first} != {log_lh_second}"
      );
    }
  }
}
