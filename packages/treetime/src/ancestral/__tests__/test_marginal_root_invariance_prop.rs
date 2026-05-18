#[cfg(test)]
mod tests {
  use crate::ancestral::__tests__::prop_generators::input::{MarginalTestInput, arb_marginal_input_no_gaps};
  use crate::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
  use proptest::prelude::*;

  proptest! {
    #![proptest_config(ProptestConfig::with_cases(30))]

    /// Felsenstein's pulley principle: for a reversible GTR model, the total log-likelihood
    /// is invariant to root placement on an unrooted tree. Dense representation.
    ///
    /// Generates a random tree+alignment+GTR, reroots at a different internal node
    /// (preserving the unrooted topology and branch lengths), and verifies that
    /// dense marginal log-likelihoods agree between the two rootings.
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

    /// Felsenstein's pulley principle: for a reversible GTR model, the total log-likelihood
    /// is invariant to root placement on an unrooted tree. Sparse representation.
    ///
    /// Same strategy as the dense test, but uses sparse marginal reconstruction
    /// (Fitch compression + marginal on variable positions only).
    ///
    /// Sparse root invariance: Fitch forward pass resolves ambiguous state
    /// sets using parent states, which differ under different rootings. This
    /// produces different compression patterns and different sets of edge
    /// mutations, causing root-dependent likelihood contributions.
    /// See kb/issues/M-ancestral-sparse-root-invariance.md.
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
    use treetime_graph::reroot::{apply_reroot_topology, remove_node_if_trivial};
    use crate::payload::ancestral::GraphAncestral;
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};

    /// Reroot a tree at a non-root internal node and return the new Newick string.
    ///
    /// Inverts edges along the path from the selected internal node to the old
    /// root, then collapses the old root (now degree-2) by merging its two edges
    /// into one with summed branch length. This preserves the unrooted topology.
    pub fn reroot_at_internal_node(newick: &str, node_idx: usize) -> Result<String, Report> {
      let mut graph: GraphAncestral = nwk_read_str(newick)?;

      let old_root_key = graph.get_exactly_one_root()?.read_arc().key();

      let internal_keys: Vec<GraphNodeKey> = graph
        .get_nodes()
        .iter()
        .filter_map(|node| {
          let node = node.read_arc();
          (!node.is_root() && !node.is_leaf()).then_some(node.key())
        })
        .sorted()
        .collect();

      if internal_keys.is_empty() {
        return Err(eyre::eyre!("No non-root internal nodes available for rerooting"));
      }

      let new_root_key = internal_keys[node_idx % internal_keys.len()];

      apply_reroot_topology(&mut graph, old_root_key, new_root_key)?;
      remove_node_if_trivial(&mut graph, old_root_key)?;

      let options = NwkWriteOptions {
        weight_significant_digits: Some(17),
        ..NwkWriteOptions::default()
      };
      nwk_write_str(&graph, &options)
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

        let graph: GraphAncestral = nwk_read_str(&rerooted)?;
        let leaves = graph.get_leaves();
        assert_eq!(leaves.len(), 4, "Must have 4 leaves: {rerooted}");

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
