#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input_small;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{
    run_dense_marginal_for_input, run_sparse_marginal_for_input,
  };
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::seq::composition::Composition;
  use proptest::prelude::*;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::prop_assert_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn test_prop_marginal_idempotency_dense(input in arb_marginal_input_small()) {
      let nwk_parsed = nwk_read(input.newick.as_bytes()).unwrap();
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let (_, recon) = run_dense_marginal_for_input(&input).unwrap();

      let (recon, log_lh_first) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths).unwrap();
      let (_, log_lh_second) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths).unwrap();

      prop_assert_abs_diff_eq!(log_lh_first.value(), log_lh_second.value(), epsilon = 1e-10);
    }

    #[test]
    fn test_prop_marginal_idempotency_sparse(input in arb_marginal_input_small()) {
      let nwk_parsed = nwk_read(input.newick.as_bytes()).unwrap();
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let (_, recon) = run_sparse_marginal_for_input(&input).unwrap();

      let (recon, log_lh_first) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths).unwrap();
      let (_, log_lh_second) = recon.marginal_update(&nwk_parsed.graph, &branch_lengths).unwrap();

      prop_assert_abs_diff_eq!(log_lh_first.value(), log_lh_second.value(), epsilon = 1e-10);
    }

    #[test]
    fn test_prop_marginal_sparse_map_composition_matches_sequence(input in arb_marginal_input_small()) {
      let (_, recon) = run_sparse_marginal_for_input(&input).unwrap();

      let compositions_match = recon.node_states.iter().all(|(key, node)| {
        let expected = Composition::with_seq(&**node.sequence, &recon.partition.alphabet);
        expected == recon.partition.obs_nodes[key].composition
      });
      prop_assert!(compositions_match);
    }
  }
}
