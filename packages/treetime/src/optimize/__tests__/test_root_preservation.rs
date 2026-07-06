#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::dispatch::run_optimize_mixed_inner;
  use crate::optimize::params::BranchOptMethod;
  use crate::partition::marginal_dense::PartitionMarginalDense;
  use crate::partition::traits::PartitionOptimizeOps;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  fn setup_dense(
    newick: &str,
    fasta: &str,
  ) -> Result<(GraphAncestral, Vec<Arc<RwLock<dyn PartitionOptimizeOps>>>), Report> {
    let alphabet = Alphabet::new(AlphabetName::Nuc)?;
    let aln = read_many_fasta_str(fasta, &alphabet)?;
    let graph: GraphAncestral = nwk_read_str(newick)?;
    let gtr = jc69(JC69Params::default())?;
    let partition = PartitionMarginalDense::new(0, gtr, alphabet, get_common_length(&aln)?);
    let partitions: Vec<Arc<RwLock<PartitionMarginalDense>>> = vec![Arc::new(RwLock::new(partition))];
    initialize_marginal(&graph, &partitions, &aln)?;
    update_marginal(&graph, &partitions)?;
    let mixed: Vec<Arc<RwLock<dyn PartitionOptimizeOps>>> = partitions
      .into_iter()
      .map(|p| -> Arc<RwLock<dyn PartitionOptimizeOps>> { p })
      .collect();
    Ok((graph, mixed))
  }

  fn root_edge_branch_lengths(graph: &GraphAncestral) -> (f64, f64) {
    let root = graph.get_exactly_one_root().unwrap();
    let children = graph.children_of(&root.read_arc());
    assert_eq!(children.len(), 2);
    let bl0 = children[0]
      .1
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap_or(0.0);
    let bl1 = children[1]
      .1
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap_or(0.0);
    (bl0, bl1)
  }

  #[rustfmt::skip]
  const BIFURCATING_TREE: &str = "((A:0.05,B:0.05)AB:0.15,(C:0.05,D:0.05)CD:0.05)root:0.0;";

  // Subtree AB has A/G at variable sites, subtree CD has T/C.
  // This ensures the root arc carries real substitution signal and is not zero-optimal.
  const ALIGNMENT: &str = indoc! {r#"
    >A
    ACGTACGTAAGGACGT
    >B
    ACGTACGTAAGGACGA
    >C
    TCGATCGATTCCACGT
    >D
    TCGATCGATTCCACGC
  "#};

  #[test]
  fn test_root_preservation_ratio_preserved_after_optimization() -> Result<(), Report> {
    let (graph, partitions) = setup_dense(BIFURCATING_TREE, ALIGNMENT)?;

    let (bl0_before, bl1_before) = root_edge_branch_lengths(&graph);
    let total_before = bl0_before + bl1_before;
    let ratio_before = bl0_before / total_before;

    run_optimize_mixed_inner(&graph, &partitions, BranchOptMethod::BrentSqrt, 0.0, true)?;

    let (bl0_after, bl1_after) = root_edge_branch_lengths(&graph);
    let total_after = bl0_after + bl1_after;
    let ratio_after = bl0_after / total_after;

    assert_abs_diff_eq!(ratio_before, ratio_after, epsilon = 1e-12);
    assert!(total_after > 0.0, "optimized total should be positive");
    Ok(())
  }

  #[test]
  fn test_root_preservation_both_edges_zero_uses_equal_split() -> Result<(), Report> {
    let tree = "((A:0.1,B:0.2)AB:0.0,(C:0.15,D:0.12)CD:0.0)root:0.0;";
    let (graph, partitions) = setup_dense(tree, ALIGNMENT)?;

    run_optimize_mixed_inner(&graph, &partitions, BranchOptMethod::BrentSqrt, 0.0, true)?;

    let (bl0, bl1) = root_edge_branch_lengths(&graph);
    let total = bl0 + bl1;
    assert!(
      total > 0.0,
      "optimizer should produce positive total from zero-start edges"
    );
    let ratio = bl0 / total;
    assert_abs_diff_eq!(ratio, 0.5, epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_root_preservation_skips_trifurcating_root() -> Result<(), Report> {
    let tree = "((A:0.05,B:0.05)AB:0.1,(C:0.05,D:0.05)CD:0.1,(E:0.05,F:0.05)EF:0.1)root:0.0;";
    let fasta = indoc! {r#"
      >A
      AAAAACGTACGTACGT
      >B
      AAAAACGTACGTACGA
      >C
      CCCCACGTACGTACGG
      >D
      CCCCACGTACGTACGC
      >E
      GGGGACGTACGTACGT
      >F
      GGGGACGTACGTACGA
    "#};
    let (graph, partitions) = setup_dense(tree, fasta)?;

    {
      let root = graph.get_exactly_one_root()?;
      let children = graph.children_of(&root.read_arc());
      assert_eq!(children.len(), 3);
    }

    // Trifurcating root: redistribution is skipped (only applies to len==2).
    // The function should complete without error.
    run_optimize_mixed_inner(&graph, &partitions, BranchOptMethod::BrentSqrt, 0.0, true)?;
    Ok(())
  }

  #[test]
  fn test_root_preservation_asymmetric_ratio() -> Result<(), Report> {
    let tree = "((A:0.05,B:0.05)AB:0.9,(C:0.05,D:0.05)CD:0.01)root:0.0;";
    let (graph, partitions) = setup_dense(tree, ALIGNMENT)?;

    let (bl0_before, bl1_before) = root_edge_branch_lengths(&graph);
    let total_before = bl0_before + bl1_before;
    let ratio_before = bl0_before / total_before;

    assert!(ratio_before > 0.9, "pre-condition: asymmetric ratio");

    run_optimize_mixed_inner(&graph, &partitions, BranchOptMethod::BrentSqrt, 0.0, true)?;

    let (bl0_after, bl1_after) = root_edge_branch_lengths(&graph);
    let total_after = bl0_after + bl1_after;
    let ratio_after = bl0_after / total_after;

    assert_abs_diff_eq!(ratio_before, ratio_after, epsilon = 1e-12);
    Ok(())
  }
}
