#[cfg(test)]
mod tests {
  use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
  use crate::commands::clock::run::run_clock;
  use crate::commands::shared::clock_filter::ClockFilterArgs;
  use crate::commands::shared::metadata::MetadataArgs;
  use crate::commands::shared::output_args::{ClockOutputSelection, ClockOutputSelectionArgs, OutputCoreArgs};
  use crate::commands::shared::reroot::RootingArgs;
  use crate::commands::shared::tree_input::TreeArgs;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::fs;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::vec_of_owned;

  const TREE: &str = "((A:0.1,B:0.2):0.1,((C:0.3,D:0.4):0.1,((E:0.2,F:0.3):0.1,(G:0.25,H:0.2):0.1):0.1));";

  const METADATA: &str =
    "strain\tdate\nA\t2016-01-01\nB\t2017-01-01\nC\t2018-01-01\nD\t2019-01-01\nE\tXXXX-XX-XX\nF\t\nG\t2019-XX-XX\n";

  #[rustfmt::skip]
  #[rstest]
  #[case::keeps_every_tip_by_default( false, vec_of_owned!["A", "B", "C", "D", "E", "F", "G", "H"])]
  #[case::removes_unknown_empty_and_missing_dates( true,  vec_of_owned!["A", "B", "C", "D", "G"])]
  #[trace]
  fn test_clock_prune_undated_tips_removes_tips_without_a_usable_date(
    #[case] prune_undated_tips: bool,
    #[case] expected_leaves: Vec<String>,
  ) -> Result<(), Report> {
    let dir = tempfile::tempdir()?;
    let tree = dir.path().join("tree.nwk");
    let metadata = dir.path().join("metadata.tsv");
    fs::write(&tree, format!("{TREE}\n"))?;
    fs::write(&metadata, METADATA)?;
    let raw = TreetimeClockArgsRaw {
      tree: TreeArgs { tree: Some(tree) },
      metadata: MetadataArgs { metadata: Some(metadata) },
      clock_filter: ClockFilterArgs {
        clock_filter: 0.0,
        prune_undated_tips,
        ..ClockFilterArgs::default()
      },
      rooting: RootingArgs {
        keep_root: true,
        ..RootingArgs::default()
      },
      output: OutputCoreArgs {
        output_all: Some(dir.path().join("output")),
        ..OutputCoreArgs::default()
      },
      output_selection: ClockOutputSelectionArgs {
        output_selection: vec![ClockOutputSelection::Nwk],
      },
      ..TreetimeClockArgsRaw::default()
    };

    let result = run_clock(&TreetimeClockArgs::try_from(raw)?, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let leaves = result
      .regression_results
      .into_iter()
      .filter(|row| row.is_leaf)
      .map(|row| row.name)
      .sorted()
      .collect_vec();
    assert_eq!(expected_leaves, leaves);
    Ok(())
  }
}
