#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::__tests__::prop_generators::input::{
    arb_discrete_input, arb_marginal_input_gap_runs, arb_marginal_input_small,
  };
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::constants::MIN_BRANCH_LENGTH_FRACTION;
  use crate::gtr::gtr::GTR;
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::discrete::partition::PartitionMarginalDiscrete;
  use crate::partition::marginal::shared::update::MarginalPasses;
  use crate::partition::storage::discrete::DiscreteStates;
  use crate::progress::NoopProgress;
  use crate::seq::alignment::node_seq_inputs;
  use crate::seq::composition::Composition;
  use ndarray::array;
  use proptest::prelude::*;

  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlphabetLike;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
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
      let (_, recon) = run_sparse_marginal(&input).unwrap();

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

    #[test]
    fn test_prop_marginal_dense_update_idempotent(input in arb_marginal_input_gap_runs()) {
      let nwk_parsed = nwk_read_str(&input.newick).unwrap();
      let names = nwk_parsed.names();
      let graph: Graph = nwk_parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let partition = PartitionMarginalDense::new(
        0,
        Alphabet::new(AlphabetName::Nuc).unwrap(),
        &graph,
        &node_seq_inputs(&graph, &names, input.alignment.clone()),
      )
      .unwrap();

      let first = partition.marginal_update(&input.gtr, &graph, &branch_lengths, &()).unwrap();
      let second = partition.marginal_update(&input.gtr, &graph, &branch_lengths, &()).unwrap();

      prop_assert_eq!(
        json_write_str(&first, JsonPretty(false)).unwrap(),
        json_write_str(&second, JsonPretty(false)).unwrap()
      );
    }

    #[test]
    fn test_prop_marginal_discrete_update_idempotent((newick, traits) in arb_discrete_input()) {
      let nwk_parsed = nwk_read_str(&newick).unwrap();
      let names = nwk_parsed.names();
      let graph: Graph = nwk_parsed.graph;
      let branch_lengths = branch_lengths_or_zero(&nwk_parsed.branch_lengths);
      let states = DiscreteStates::from_values(["a", "b", "c"].into_iter(), "?");
      let gtr = GTR::builder().n_states(3).mu(1.0).pi(array![0.2, 0.3, 0.5]).build().unwrap();
      let partition = PartitionMarginalDiscrete::new(
        states,
        &graph,
        &traits,
        &names,
        MIN_BRANCH_LENGTH_FRACTION,
        false,
        &NoopProgress,
      )
      .unwrap();

      let first = partition.marginal_update(&gtr, &graph, &branch_lengths, &()).unwrap();
      let second = partition.marginal_update(&gtr, &graph, &branch_lengths, &()).unwrap();

      prop_assert_eq!(
        json_write_str(&first, JsonPretty(false)).unwrap(),
        json_write_str(&second, JsonPretty(false)).unwrap()
      );
    }
  }
}
