#[cfg(test)]
mod tests {
  use crate::commands::ancestral::marginal::update_marginal;
  use crate::commands::optimize::partition_ops::PartitionOptimizeOps;
  use crate::representation::partition::traits::ExactStateCache;
  use crate::representation::payload::ancestral::GraphAncestral;
  use approx::assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

  use super::super::test_dense_sparse_equivalence_support::tests::{
    NUC_ALPHABET, TREE_NEWICK, setup_dense_only, setup_sparse_only,
  };

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence() -> Result<(), Report> {
    let aln = super::super::test_dense_sparse_equivalence_support::tests::gap_free_alignment()?;

    // Initialize dense
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;
    let log_lh_dense = update_marginal(&graph_dense, &dense_partitions)?;

    // Initialize sparse
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;
    let log_lh_sparse = update_marginal(&graph_sparse, &sparse_partitions)?;

    // Initial log-LH should be equivalent (before any optimization)
    assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence_with_mutations() -> Result<(), Report> {
    // Alignment with more mutations
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      AAAAAAAAAAAAAAAA
      >B
      CCCCCCCCCCCCCCCC
      >C
      GGGGGGGGGGGGGGGG
      >D
      TTTTTTTTTTTTTTTT
    "#},
      &*NUC_ALPHABET,
    )?;

    // Initialize dense
    let graph_dense: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;
    let log_lh_dense = update_marginal(&graph_dense, &dense_partitions)?;

    // Initialize sparse
    let graph_sparse: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;
    let log_lh_sparse = update_marginal(&graph_sparse, &sparse_partitions)?;

    // Initial log-LH should be equivalent
    assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }
  #[test]
  fn test_dense_sparse_ambiguous_r_edge_contribution_equivalence() -> Result<(), Report> {
    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      GCCC
      >B
      RCCC
      >C
      ACCC
      >D
      ACCC
    "#},
      &*NUC_ALPHABET,
    )?;

    let graph_dense: GraphAncestral = nwk_read_str("((A:0.05,B:0.05)AB:0.05,(C:0.05,D:0.05)CD:0.05)root:0.01;")?;
    let dense_partitions = setup_dense_only(&graph_dense, &aln)?;
    update_marginal(&graph_dense, &dense_partitions)?;

    let graph_sparse: GraphAncestral = nwk_read_str("((A:0.05,B:0.05)AB:0.05,(C:0.05,D:0.05)CD:0.05)root:0.01;")?;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &aln)?;

    let dense_edge_key = graph_dense
      .get_edges()
      .iter()
      .find_map(|edge_ref| {
        let edge = edge_ref.read_arc();
        let child_name = graph_dense
          .get_node(edge.target())
          .expect("child exists")
          .read_arc()
          .payload()
          .read_arc()
          .name
          .clone();
        (child_name.as_deref() == Some("B")).then_some(edge.key())
      })
      .expect("edge to B exists");
    let sparse_edge_key = graph_sparse
      .get_edges()
      .iter()
      .find_map(|edge_ref| {
        let edge = edge_ref.read_arc();
        let child_name = graph_sparse
          .get_node(edge.target())
          .expect("child exists")
          .read_arc()
          .payload()
          .read_arc()
          .name
          .clone();
        (child_name.as_deref() == Some("B")).then_some(edge.key())
      })
      .expect("edge to B exists");

    let mut dense_cache = ExactStateCache::new();
    let dense_contribution =
      dense_partitions[0]
        .read_arc()
        .create_edge_contribution(&graph_dense, dense_edge_key, &mut dense_cache)?;
    let mut sparse_cache = ExactStateCache::new();
    let sparse_contribution =
      sparse_partitions[0]
        .read_arc()
        .create_edge_contribution(&graph_sparse, sparse_edge_key, &mut sparse_cache)?;

    for branch_length in [0.01, 0.05, 0.1] {
      let dense_metrics = dense_contribution.evaluate(branch_length);
      let sparse_metrics = sparse_contribution.evaluate(branch_length);

      assert_ulps_eq!(dense_metrics.log_lh, sparse_metrics.log_lh, max_ulps = 128);
      assert_ulps_eq!(dense_metrics.derivative, sparse_metrics.derivative, max_ulps = 128);
      assert_ulps_eq!(
        dense_metrics.second_derivative,
        sparse_metrics.second_derivative,
        max_ulps = 128
      );
    }

    Ok(())
  }
}
