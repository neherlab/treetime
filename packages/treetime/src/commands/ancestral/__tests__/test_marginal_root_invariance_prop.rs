#[cfg(test)]
mod tests {
  use crate::commands::ancestral::__tests__::prop_generators::input::{MarginalTestInput, arb_marginal_input_no_gaps};
  use crate::commands::ancestral::__tests__::prop_marginal_support::tests::{run_dense_marginal, run_sparse_marginal};
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
    /// Tolerance is looser than dense (1e-1 vs 1e-6) because Fitch forward pass
    /// resolves ambiguous state sets using parent states, which differ under
    /// different rootings. This produces different compression patterns and
    /// different sets of edge mutations, causing root-dependent likelihood
    /// contributions beyond the matrix exponential path difference.
    /// Measured max diff: ~3e-4 typical, ~1e-2 with rare proptest seeds;
    /// 1e-1 is tightest 1e-N that accommodates the worst observed case.
    #[test]
    fn test_prop_marginal_sparse_log_lh_root_invariance(
      input in arb_marginal_input_no_gaps(4, 10),
      node_idx in 0_usize..100,
    ) {
      let rerooted_newick = helpers::reroot_at_internal_node(&input.newick, node_idx).unwrap();

      let lh1 = run_sparse_marginal(&input).unwrap().0;
      let input2 = MarginalTestInput { newick: rerooted_newick, ..input };
      let lh2 = run_sparse_marginal(&input2).unwrap().0;

      let diff = (lh1 - lh2).abs();
      prop_assert!(diff < 1e-1,
        "Sparse root invariance violated: lh1={lh1}, lh2={lh2}, diff={diff}");
    }
  }

  mod helpers {
    use crate::representation::payload::ancestral::GraphAncestral;
    use eyre::Report;
    use itertools::Itertools;
    use treetime_graph::edge::{GraphEdgeKey, HasBranchLength, invert_edge};
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

      let path = graph.path_from_node_to_node(new_root_key, old_root_key)?;
      for (_, edge) in &path {
        if let Some(edge) = edge {
          invert_edge(&mut graph, edge);
        }
      }
      drop(path);

      // Bypass the degree-2 old root, then remove it and the dangling edge
      // so the graph is in a consistent state (no ghost nodes or orphaned edges).
      let inbound_key = bypass_degree2_node(&graph, old_root_key)?;
      graph.remove_edge(inbound_key)?;
      graph.remove_node(old_root_key)?;

      graph.build()?;

      let options = NwkWriteOptions {
        weight_significant_digits: Some(17),
        ..NwkWriteOptions::default()
      };
      nwk_write_str(&graph, &options)
    }

    /// Bypass a degree-2 node (1 inbound, 1 outbound) by redirecting its
    /// outbound edge to connect parent directly to child, with summed branch length.
    ///
    /// Returns the inbound edge key so the caller can remove it (along with the
    /// now-disconnected node) to leave the graph consistent.
    fn bypass_degree2_node(graph: &GraphAncestral, node_key: GraphNodeKey) -> Result<GraphEdgeKey, Report> {
      let node = graph
        .get_node(node_key)
        .ok_or_else(|| eyre::eyre!("Node {node_key} not found"))?;
      let node_ref = node.read_arc();

      if node_ref.inbound().len() != 1 || node_ref.outbound().len() != 1 {
        return Err(eyre::eyre!(
          "Node {node_key} is not degree-2: in={}, out={}",
          node_ref.inbound().len(),
          node_ref.outbound().len()
        ));
      }

      let inbound_key = node_ref.inbound()[0];
      let outbound_key = node_ref.outbound()[0];
      drop(node_ref);

      let inbound_edge = graph
        .get_edge(inbound_key)
        .ok_or_else(|| eyre::eyre!("Inbound edge {inbound_key} not found"))?;
      let outbound_edge = graph
        .get_edge(outbound_key)
        .ok_or_else(|| eyre::eyre!("Outbound edge {outbound_key} not found"))?;

      let parent_key = inbound_edge.read_arc().source();
      let inbound_bl = inbound_edge
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap_or(0.0);
      let outbound_bl = outbound_edge
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap_or(0.0);

      // Redirect outbound edge: old_root->child becomes parent->child
      outbound_edge.write_arc().set_source(parent_key);
      outbound_edge
        .write_arc()
        .payload()
        .write_arc()
        .set_branch_length(Some(inbound_bl + outbound_bl));

      // Update parent: replace inbound_key with outbound_key in outbound list
      let parent = graph
        .get_node(parent_key)
        .ok_or_else(|| eyre::eyre!("Parent node {parent_key} not found"))?;
      {
        let mut parent_ref = parent.write_arc();
        parent_ref.outbound_mut().retain(|&k| k != inbound_key);
        parent_ref.outbound_mut().push(outbound_key);
      }

      // Clear old root's outbound so remove_node won't find the redirected edge
      node.write_arc().outbound_mut().clear();

      Ok(inbound_key)
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
