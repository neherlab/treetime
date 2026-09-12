#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input_small;
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::ancestral::marginal::{ancestral_reconstruction, profile_branch_lengths};
  use crate::ancestral::pipeline::SparseReconstruction;
  use crate::ancestral::sample::SampleMode;
  use crate::seq::composition::Composition;
  use proptest::prelude::*;
  use rand::SeedableRng;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_primitives::AlphabetLike;
  use treetime_utils::prop_assert_abs_diff_eq;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(50))]

    /// Property test: `marginal_update` is a fixed-point operation on dense partitions
    /// across randomly generated trees, alignments, and GTR parameters.
    ///
    /// Felsenstein's two-pass pruning algorithm - backward (postorder, leaves to root)
    /// followed by forward (preorder, root to leaves) - computes exact marginal
    /// posterior state probabilities at every node. On tree-structured models this is
    /// equivalent to the sum-product algorithm (belief propagation), which converges
    /// in a single message-passing schedule. The log-likelihood, computed from root
    /// partial likelihoods after the backward pass, is therefore a deterministic
    /// function of the input data and model parameters alone.
    ///
    /// Re-running `marginal_update` on partition state that already contains computed
    /// profiles must reproduce the same log-likelihood. Any deviation indicates state
    /// corruption during a pass, accumulation into node profiles without clearing
    /// previous values, or the algorithm incorrectly depending on intermediate state
    /// from a prior run.
    ///
    /// Random inputs exercise diverse tree topologies, branch lengths, sequence
    /// compositions (including ambiguity codes and gaps), equilibrium frequencies,
    /// and substitution rate matrices to ensure the invariant holds regardless of
    /// parameterization.
    ///
    /// Invariant: let `(L1, P') = marginal_update(P)` and `(L2, _) = marginal_update(P')`,
    /// then `L1 == L2`.
    ///
    /// Companion example test: `test_marginal_idempotency_example_dense`.
    #[test]
    fn test_prop_marginal_idempotency_dense(input in arb_marginal_input_small()) {
      let NwkParse { graph, branch_lengths, .. } = nwk_read_str(&input.newick).unwrap();
      let graph: Graph = graph;
      let (_, mut recon) = run_dense_marginal(&input).unwrap();

      let log_lh_first = recon.run_marginal_update(&graph, &profile_branch_lengths(&branch_lengths))
        .unwrap()
        .value();
      let log_lh_second = recon.run_marginal_update(&graph, &profile_branch_lengths(&branch_lengths))
        .unwrap()
        .value();

      prop_assert_abs_diff_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);
    }

    /// Property test: `marginal_update` is a fixed-point operation on sparse partitions
    /// across randomly generated trees, alignments, and GTR parameters.
    ///
    /// The sparse representation compresses sequences via Fitch parsimony: only
    /// variable (mutated) positions store full probability vectors, while invariant
    /// positions are grouped and handled with a separate accumulation path. This
    /// test verifies that the compression-aware backward and forward passes are
    /// idempotent, exercising the split between per-position and grouped-invariant
    /// likelihood contributions that the dense variant does not have.
    ///
    /// The same sum-product / belief propagation exactness guarantee applies: on a
    /// tree, a single two-pass schedule produces exact marginals, so re-running the
    /// algorithm on already-computed partition state must reproduce the same
    /// log-likelihood.
    ///
    /// Invariant: let `(L1, P') = marginal_update(P)` and `(L2, _) = marginal_update(P')`,
    /// then `L1 == L2`.
    ///
    /// Companion example test: `test_marginal_idempotency_example_sparse`.
    #[test]
    fn test_prop_marginal_idempotency_sparse(input in arb_marginal_input_small()) {
      let NwkParse { graph, branch_lengths, .. } = nwk_read_str(&input.newick).unwrap();
      let graph: Graph = graph;
      let (_, mut recon) = run_sparse_marginal(&input).unwrap();

      let log_lh_first = recon.run_marginal_update(&graph, &profile_branch_lengths(&branch_lengths))
        .unwrap()
        .value();
      let log_lh_second = recon.run_marginal_update(&graph, &profile_branch_lengths(&branch_lengths))
        .unwrap()
        .value();

      prop_assert_abs_diff_eq!(log_lh_first, log_lh_second, epsilon = 1e-10);
    }

    /// Each sparse node's cached composition describes its stored parsimony sequence, and MAP
    /// reconstruction leaves both untouched. `combine_messages` and `infer_gtr` both depend on this
    /// pairing, so an output pass must not rewrite either.
    #[test]
    fn test_prop_marginal_sparse_map_composition_matches_sequence(input in arb_marginal_input_small()) {
      let NwkParse { graph, branch_lengths, .. } = nwk_read_str(&input.newick).unwrap();
      let graph: Graph = graph;
      let (_, mut recon) = run_sparse_marginal(&input).unwrap();
      let mut rng = rand::rngs::StdRng::seed_from_u64(0);

      {
        let SparseReconstruction {
          partition,
          node_states,
          forward,
          ..
        } = &mut recon;
        ancestral_reconstruction(
          &graph,
          |node| partition.reconstruct_node_sequence(node_states, forward, node, true, false, SampleMode::Argmax, &mut rng),
          |_, _| Ok(()),
        )
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
