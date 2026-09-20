#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AlignmentRecord;

  use super::super::test_dense_sparse_equivalence_support::tests::{
    NUC_ALPHABET, TREE_NEWICK, setup_dense_only, setup_sparse_only,
  };

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence() -> Result<(), Report> {
    let aln = super::super::test_dense_sparse_equivalence_support::tests::gap_free_alignment()?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let branch_lengths_dense = nwk_parsed.branch_lengths;
    let dense_partitions = setup_dense_only(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let (dense_partitions, log_lh_dense) = marginal_update_dense(
      &graph_dense,
      &branch_lengths_or_zero(&branch_lengths_dense),
      dense_partitions,
    )?;
    let log_lh_dense = log_lh_dense.value();

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let branch_lengths_sparse = nwk_parsed.branch_lengths;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;
    let (sparse_partitions, log_lh_sparse) = marginal_update_sparse(
      &graph_sparse,
      &branch_lengths_or_zero(&branch_lengths_sparse),
      sparse_partitions,
    )?;
    let log_lh_sparse = log_lh_sparse.value();

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence_with_mutations() -> Result<(), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(
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
    )?
    .into_iter()
    .map(AlignmentRecord::from)
    .collect();

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_dense_names = nwk_parsed.names();
    let graph_dense = nwk_parsed.graph;
    let branch_lengths_dense = nwk_parsed.branch_lengths;
    let dense_partitions = setup_dense_only(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let (dense_partitions, log_lh_dense) = marginal_update_dense(
      &graph_dense,
      &branch_lengths_or_zero(&branch_lengths_dense),
      dense_partitions,
    )?;
    let log_lh_dense = log_lh_dense.value();

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_sparse_names = nwk_parsed.names();
    let graph_sparse = nwk_parsed.graph;
    let branch_lengths_sparse = nwk_parsed.branch_lengths;
    let sparse_partitions = setup_sparse_only(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;
    let (sparse_partitions, log_lh_sparse) = marginal_update_sparse(
      &graph_sparse,
      &branch_lengths_or_zero(&branch_lengths_sparse),
      sparse_partitions,
    )?;
    let log_lh_sparse = log_lh_sparse.value();

    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }
}
