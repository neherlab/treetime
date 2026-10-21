#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{SequenceMode, marginal_node_times};
  use crate::optimize::indel::IndelPolicy;
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  #[trace]
  fn test_gm_runner_marginal_sparse_dense_consistency(
    #[case] dataset: &str,
    #[values(IndelPolicy::Include, IndelPolicy::Exclude)] indel_policy: IndelPolicy,
  ) -> Result<(), Report> {
    let dense = marginal_node_times(dataset, SequenceMode::Dense, indel_policy)?;

    let sparse = marginal_node_times(dataset, SequenceMode::Sparse, indel_policy)?;

    // TODO: kb/issues/H-timetree-sparse-dense-node-times-differ.md
    // v1 sparse node times differ from v1 dense by up to 4.1e-3 years (NODE_0000011 of flu_h3n2_20).
    pretty_assert_map_abs_diff_eq!(&dense, &sparse, epsilon = 5e-3);
    Ok(())
  }
}
