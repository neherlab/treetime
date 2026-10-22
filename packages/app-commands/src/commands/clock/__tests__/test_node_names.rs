#[cfg(test)]
mod tests {
  use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
  use crate::commands::clock::run::run_clock;
  use crate::commands::shared::metadata::MetadataArgs;
  use crate::commands::shared::output_args::ClockOutputSelectionArgs;
  use crate::commands::shared::output_args::{ClockOutputSelection, OutputCoreArgs};
  use crate::commands::shared::tree_input::TreeArgs;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use std::path::PathBuf;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::nwk::nwk_read_file;

  #[test]
  fn test_clock_names_the_new_root_the_same_in_the_tree_and_the_regression_rows() -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/flu/h3n2/20");
    let raw = TreetimeClockArgsRaw {
      tree: TreeArgs {
        tree: Some(data.join("tree.nwk")),
      },
      metadata: MetadataArgs {
        metadata: Some(data.join("metadata.tsv")),
      },
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..OutputCoreArgs::default()
      },
      output_selection: ClockOutputSelectionArgs {
        output_selection: vec![ClockOutputSelection::Nwk],
      },
      ..TreetimeClockArgsRaw::default()
    };

    let result = run_clock(
      &TreetimeClockArgs::try_from(raw)?,
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    let tree = nwk_read_file(output.path().join("clock.nwk"))?;
    let node_names = tree
      .graph
      .get_nodes()
      .map(|node| tree.names[node.key()].to_owned())
      .sorted()
      .collect_vec();
    let row_names = result
      .regression_results
      .into_iter()
      .map(|row| row.name)
      .sorted()
      .collect_vec();
    assert_eq!(BTreeSet::new(), tree.unnamed_nodes);
    assert_eq!(node_names, row_names);
    Ok(())
  }
}
