#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed, run_optimize_mixed_with_indel_rate};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::test_utils::dense_partition_with_constant_leaves;
  use eyre::Report;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;

  #[test]
  fn test_optimize_zero_sequence_length_run_optimize_error() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;".as_slice())?;
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let total_length = dense_partition_with_constant_leaves(&graph, Alphabet::new(AlphabetName::Nuc)?, 0)?.length;
    let result = run_optimize_mixed(
      &graph,
      total_length,
      &BTreeMap::new(),
      &BTreeMap::new(),
      BranchOptMethod::Newton,
      &mut branch_lengths,
    );
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot optimize branch lengths"
    );
    Ok(())
  }

  #[test]
  fn test_optimize_zero_sequence_length_initial_guess_error() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;".as_slice())?;
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let total_length = dense_partition_with_constant_leaves(&graph, Alphabet::new(AlphabetName::Nuc)?, 0)?.length;
    let result = initial_guess_mixed(
      &graph,
      total_length,
      &BTreeMap::new(),
      &BTreeMap::new(),
      &BTreeMap::new(),
      ExistingBranchLengths::Overwrite,
      false,
      &mut branch_lengths,
    );
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot compute initial guess"
    );
    Ok(())
  }

  #[test]
  fn test_optimize_zero_sequence_length_run_optimize_with_fixed_rate_error() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,C:0.2)root:0.01;".as_slice())?;
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let total_length = dense_partition_with_constant_leaves(&graph, Alphabet::new(AlphabetName::Nuc)?, 0)?.length;
    let result = run_optimize_mixed_with_indel_rate(
      &graph,
      total_length,
      &BTreeMap::new(),
      &BTreeMap::new(),
      BranchOptMethod::Newton,
      1.0,
      &mut branch_lengths,
    );
    assert_error!(
      result,
      "Total sequence length across all partitions is zero; cannot optimize branch lengths"
    );
    Ok(())
  }
}
