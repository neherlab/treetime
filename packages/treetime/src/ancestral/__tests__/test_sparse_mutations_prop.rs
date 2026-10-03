#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::arb_marginal_input;
  use crate::ancestral::__tests__::prop_marginal_support::tests::run_sparse_marginal;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::partition::marginal::sparse::mutations::sparse_edge_mutations;
  use crate::seq::mutation::{MutationTrack, stream_sequence_mutations};
  use proptest::prelude::*;
  use treetime_io::nwk::nwk_read_str;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(300))]

    #[test]
    fn test_prop_sparse_edge_mutations_match_full_sequence_comparison(
      input in arb_marginal_input(),
      impute in any::<bool>(),
    ) {
      let graph = nwk_read_str(&input.newick).unwrap().graph;
      let (_, sparse) = run_sparse_marginal(&input).unwrap();
      let track = MutationTrack::Nucleotide;

      let actual =
        sparse_edge_mutations(&sparse.partition, &graph, &sparse.node_states, &sparse.edges.forward, impute, &track)
          .unwrap();

      let reconstruction = MarginalReconstruction::Sparse(sparse);
      let expected = stream_sequence_mutations(
        &graph,
        reconstruction.alphabet(),
        &track,
        false,
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
