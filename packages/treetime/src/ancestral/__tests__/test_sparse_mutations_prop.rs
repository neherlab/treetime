#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input;
  use crate::ancestral::__tests__::prop_marginal_support::tests::run_sparse_marginal;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::partition::marginal::sparse::mutations::sparse_edge_mutations;
  use crate::seq::mutation::{MutationTrack, stream_sequence_mutations};
  use proptest::prelude::*;
  use treetime_io::nwk::nwk_read;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]

    #[test]
    fn test_prop_sparse_edge_mutations_match_full_sequence_comparison(
      input in arb_marginal_input(),
      impute in any::<bool>(),
      without_edge_fitch_subs in any::<bool>(),
      report_unknown in any::<bool>(),
    ) {
      let graph = nwk_read(input.newick.as_bytes()).unwrap().graph;
      let (_, mut sparse) = run_sparse_marginal(&input).unwrap();
      if without_edge_fitch_subs {
        sparse.partition.obs_edges.values_mut().for_each(|edge| edge.set_fitch_subs(vec![]));
      }
      let track = MutationTrack::Nucleotide;

      let actual =
        sparse_edge_mutations(
          &sparse.partition,
          &graph,
          &sparse.node_states,
          &sparse.edges.forward,
          impute,
          report_unknown,
          &track,
        )
          .unwrap();

      let reconstruction = MarginalReconstruction::Sparse(sparse);
      let expected = stream_sequence_mutations(
        &graph,
        reconstruction.alphabet(),
        &track,
        false,
        report_unknown,
        |node_key| reconstruction.node_sequence(&graph, impute, node_key),
        |edge_key| reconstruction.edge_indels(edge_key),
        None,
      )
      .unwrap()
      .edge_mutations;

      prop_assert_eq!(expected, actual);
    }
  }
}
