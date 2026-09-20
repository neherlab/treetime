#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::{MarginalTestInput, arb_marginal_input_no_gaps};
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use proptest::prelude::*;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    #[test]
    fn test_prop_marginal_dense_log_lh_root_invariance(
      input in arb_marginal_input_no_gaps(4, 10),
      node_idx in 0_usize..100,
    ) {
      let rerooted_newick = helpers::reroot_at_internal_node(&input.newick, node_idx).unwrap();

      let lh1 = run_dense_marginal(&input).unwrap().0;
      let input2 = MarginalTestInput { newick: rerooted_newick, ..input };
      let lh2 = run_dense_marginal(&input2).unwrap().0;

      let diff = (lh1 - lh2).abs();
      prop_assert!(diff < 1e-6,
        "Dense root invariance violated: lh1={lh1}, lh2={lh2}, diff={diff}");
    }

    #[test]
    #[ignore = "sparse root invariance violation: max ~1e-2 (kb/issues/M-ancestral-sparse-root-invariance.md)"]
    fn test_prop_marginal_sparse_log_lh_root_invariance(
      input in arb_marginal_input_no_gaps(4, 10),
      node_idx in 0_usize..100,
    ) {
      let rerooted_newick = helpers::reroot_at_internal_node(&input.newick, node_idx).unwrap();

      let lh1 = run_sparse_marginal(&input).unwrap().0;
      let input2 = MarginalTestInput { newick: rerooted_newick, ..input };
      let lh2 = run_sparse_marginal(&input2).unwrap().0;

      let diff = (lh1 - lh2).abs();
      prop_assert!(diff < 1e-6,
        "Sparse root invariance violated: lh1={lh1}, lh2={lh2}, diff={diff}");
    }
  }

  mod helpers {
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::reroot::{
      apply_reroot_topology, record_merge, remove_node_if_trivial, trivial_node_branch_lengths,
    };
    use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

    pub fn reroot_at_internal_node(newick: &str, node_idx: usize) -> Result<String, Report> {
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;
      let mut graph: Graph = graph;

      let old_root_key = graph.get_exactly_one_root()?.key();

      let internal_keys: Vec<GraphNodeKey> = graph
        .get_nodes()
        .filter_map(|node| (!node.is_root() && !node.is_leaf()).then_some(node.key()))
        .sorted()
        .collect();

      if internal_keys.is_empty() {
        return Err(eyre::eyre!("No non-root internal nodes available for rerooting"));
      }

      let new_root_key = internal_keys[node_idx % internal_keys.len()];

      apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;
      let (old_root_parent, old_root_child) = trivial_node_branch_lengths(&graph, old_root_key, &branch_lengths);
      if let Some(info) = remove_node_if_trivial(&mut graph, old_root_key, old_root_parent, old_root_child)? {
        record_merge(&mut branch_lengths, &info);
      }

      let options = NwkWriteOptions {
        weight_significant_digits: Some(17),
        ..NwkWriteOptions::default()
      };
      nwk_write_str(&graph, &names, &branch_lengths, &options)
    }

    #[cfg(test)]
    mod tests {
      use super::*;

      #[test]
      fn test_reroot_at_internal_node_preserves_topology() -> Result<(), Report> {
        let newick = "((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root:0.001;";
        let rerooted = reroot_at_internal_node(newick, 0)?;

        assert!(rerooted.ends_with(';'), "Must end with semicolon: {rerooted}");

        for taxon in &["A", "B", "C", "D"] {
          assert!(rerooted.contains(taxon), "Missing taxon {taxon} in {rerooted}");
        }

        let nwk_parsed = nwk_read_str(&rerooted)?;
        let names = nwk_parsed.names();
        let graph = nwk_parsed.graph;
        let branch_lengths = nwk_parsed.branch_lengths;

        let graph: Graph = graph;
        let leaves = graph.get_leaves().collect::<Vec<_>>();
        assert_eq!(4, leaves.len(), "Must have 4 leaves: {rerooted}");

        Ok(())
      }

      #[test]
      fn test_reroot_at_internal_node_different_rootings_differ() -> Result<(), Report> {
        let newick = "((A:0.1,B:0.2)AB:0.3,(C:0.15,D:0.25)CD:0.4)root:0.001;";
        let rerooted0 = reroot_at_internal_node(newick, 0)?;
        let rerooted1 = reroot_at_internal_node(newick, 1)?;

        assert_ne!(
          rerooted0, rerooted1,
          "Different indices should produce different rootings"
        );

        Ok(())
      }
    }
  }
}
