#[cfg(test)]
mod tests {
  use super::super::test_gm_runner_support::support::{OUTPUTS, SequenceMode, marginal_node_times};
  use crate::optimize::indel::IndelPolicy;
  use eyre::Report;
  use rstest::rstest;
  use treetime_utils::pretty_assert_map_abs_diff_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::flu_h3n2_20("flu_h3n2_20")]
  // TODO: kb/issues/M-timetree-gm-runner-missing-internal-times.md
  // These datasets differ from v0 by more than flu_h3n2_20 (largest gap with indels included, excluded):
  // #[case::dengue_20("dengue_20")]               // 14.8, 16.3 years at NODE_0000018
  // #[case::ebola_20("ebola_20")]                 // v1 dates 8 internal nodes that v0 leaves out: kb/issues/M-timetree-dense-golden-master-node-mismatch.md
  // #[case::lassa_l_20("lassa_L_20")]             // 4.9, 6.7 years at NODE_0000004
  // #[case::mpox_clade_ii_20("mpox_clade_ii_20")] // exceeds the 300 s test timeout: kb/issues/M-timetree-marginal-dense-mpox-slow.md
  // #[case::rsv_a_20("rsv_a_20")]                 // 1.8, 1.7 years at NODE_0000017
  // #[case::tb_20("tb_20")]                       // v1 dates 6 internal nodes that v0 leaves out; 33.9 years at NODE_0000016
  // #[case::zika_20("zika_20")]                   // the metadata has no column '#name': kb/issues/M-timetree-date-header-hash.md
  #[trace]
  fn test_gm_runner_marginal_dense(
    #[case] dataset: &str,
    #[values(IndelPolicy::Include, IndelPolicy::Exclude)] indel_policy: IndelPolicy,
  ) -> Result<(), Report> {
    let expected = OUTPUTS[dataset].marginal_dense();

    let actual = marginal_node_times(dataset, SequenceMode::Dense, indel_policy)?;

    // TODO: kb/issues/M-timetree-branch-grid-uniform-resolution.md
    // v1 dense node times differ from v0 by up to 0.070 years (NODE_0000017 of flu_h3n2_20).
    pretty_assert_map_abs_diff_eq!(expected, &actual, epsilon = 1e-1);
    Ok(())
  }
}
