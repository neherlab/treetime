#[cfg(test)]
mod tests {
  use crate::commands::ancestral::marginal::update_marginal;
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
}
