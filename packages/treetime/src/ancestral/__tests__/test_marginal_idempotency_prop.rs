#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input_small;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::ancestral::sample::SampleMode;
  use crate::ancestral::tip_states::TipStates;
  use crate::seq::composition::Composition;
  use proptest::prelude::*;
  use rand::SeedableRng;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlphabetLike;
  use treetime_utils::prop_assert_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    #[test]
    fn test_prop_marginal_idempotency_dense(input in arb_marginal_input_small()) {
      let nwk_parsed = nwk_read_str(&input.newick).unwrap();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let (_, recon) = run_dense_marginal(&input).unwrap();

      let (recon, log_lh_first) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths)).unwrap();
      let log_lh_first = log_lh_first.value();
      let (recon, log_lh_second) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths)).unwrap();
      let log_lh_second = log_lh_second.value();

      prop_assert_abs_diff_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);
    }

    #[test]
    fn test_prop_marginal_idempotency_sparse(input in arb_marginal_input_small()) {
      let nwk_parsed = nwk_read_str(&input.newick).unwrap();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let (_, recon) = run_sparse_marginal(&input).unwrap();

      let (recon, log_lh_first) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths)).unwrap();
      let log_lh_first = log_lh_first.value();
      let (recon, log_lh_second) = recon.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths)).unwrap();
      let log_lh_second = log_lh_second.value();

      prop_assert_abs_diff_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);
    }

    #[test]
    fn test_prop_marginal_sparse_map_composition_matches_sequence(input in arb_marginal_input_small()) {
      let nwk_parsed = nwk_read_str(&input.newick).unwrap();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let (_, mut recon) = run_sparse_marginal(&input).unwrap();
      let mut rng = rand::rngs::StdRng::seed_from_u64(0);

      {
        let SparseReconstruction {
          partition,
          node_states,
          edges,
          ..

        } = &mut recon;
        ancestral_reconstruction(&graph, |node| {
          partition
            .reconstruct_node_sequence(node_states, &edges.forward, node, TipStates { include_leaves: true, impute: false }, SampleMode::Argmax, &mut rng)
            .map(|seq| seq.is_some())
        })
        .unwrap();
      }

      let compositions_match = recon.node_states.iter().all(|(key, node)| {
        let expected = Composition::with_seq(
          &node.sequence,
          recon.partition.alphabet.chars(),
          recon.partition.alphabet.gap(),
        );
        expected == recon.partition.obs_nodes[key].composition
      });
      prop_assert!(compositions_match);
    }
  }
}
