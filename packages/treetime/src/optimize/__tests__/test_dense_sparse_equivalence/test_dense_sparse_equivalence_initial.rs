#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::profile_branch_lengths;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use indoc::indoc;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::{NwkParse, nwk_read_str};

  use super::super::test_dense_sparse_equivalence_support::tests::{
    NUC_ALPHABET, TREE_NEWICK, setup_dense_only, setup_sparse_only,
  };

  #[test]
  fn test_dense_sparse_initial_log_lh_equivalence() -> Result<(), Report> {
    let aln = super::super::test_dense_sparse_equivalence_support::tests::gap_free_alignment()?;

    // Initialize dense
    let NwkParse {
      graph: graph_dense,
      names: graph_dense_names,
      branch_lengths: branch_lengths_dense,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut dense_partitions = setup_dense_only(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let log_lh_dense = marginal_update(
      &graph_dense,
      &profile_branch_lengths(&branch_lengths_dense),
      &mut dense_partitions,
    )?
    .value();

    // Initialize sparse
    let NwkParse {
      graph: graph_sparse,
      names: graph_sparse_names,
      branch_lengths: branch_lengths_sparse,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut sparse_partitions = setup_sparse_only(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;
    let log_lh_sparse = marginal_update(
      &graph_sparse,
      &profile_branch_lengths(&branch_lengths_sparse),
      &mut sparse_partitions,
    )?
    .value();

    // Initial log-LH should be equivalent (before any optimization)
    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

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
    let NwkParse {
      graph: graph_dense,
      names: graph_dense_names,
      branch_lengths: branch_lengths_dense,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut dense_partitions = setup_dense_only(&graph_dense, &graph_dense_names, &aln, &branch_lengths_dense)?;
    let log_lh_dense = marginal_update(
      &graph_dense,
      &profile_branch_lengths(&branch_lengths_dense),
      &mut dense_partitions,
    )?
    .value();

    // Initialize sparse
    let NwkParse {
      graph: graph_sparse,
      names: graph_sparse_names,
      branch_lengths: branch_lengths_sparse,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let mut sparse_partitions = setup_sparse_only(&graph_sparse, &graph_sparse_names, &aln, &branch_lengths_sparse)?;
    let log_lh_sparse = marginal_update(
      &graph_sparse,
      &profile_branch_lengths(&branch_lengths_sparse),
      &mut sparse_partitions,
    )?
    .value();

    // Initial log-LH should be equivalent
    pretty_assert_ulps_eq!(log_lh_dense, log_lh_sparse, max_ulps = 100);

    Ok(())
  }
}
