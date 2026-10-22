#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output_args::TimetreeOutputSelectionArgs;
  use crate::commands::shared::output_args::{OutputCoreArgs, TimetreeOutputSelection};
  use crate::commands::shared::seed::SeedArgs;
  use crate::commands::shared::tree_input::TreeArgs;
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use app_output::output_plan::OutputSelection;
  use app_output::table_output::table_read_file;
  use deser::Deserialize;
  use eyre::Report;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeSet;
  use std::path::PathBuf;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_io::fasta::fasta_read_file;
  use treetime_io::nwk::nwk_read_file;

  #[rstest]
  #[case::rerooted(true)]
  #[case::rerooted_and_resolved(false)]
  #[trace]
  fn test_timetree_names_every_node_the_same_in_every_output(#[case] keep_polytomies: bool) -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let data = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
    let raw = TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![data.join("aln.fasta.xz")],
      },
      tree: TreeArgs {
        tree: Some(data.join("tree.nwk")),
      },
      metadata: Some(data.join("metadata.tsv")),
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..OutputCoreArgs::default()
      },
      output_selection: TimetreeOutputSelectionArgs {
        output_selection: vec![
          TimetreeOutputSelection::Nwk,
          TimetreeOutputSelection::ConfidenceTsv,
          TimetreeOutputSelection::Sequences,
        ],
      },
      confidence: true,
      keep_polytomies,
      seed_args: SeedArgs { seed: Some(7) },
      ..TreetimeTimetreeArgsRaw::default()
    };
    run_timetree_estimation(
      &TreetimeTimetreeArgs::try_from(raw)?,
      &NoopCancel,
      &NoopProgress,
      &NoopProgress,
    )?;

    let tree = nwk_read_file(output.path().join("timetree.nwk"))?;
    let node_names = tree
      .graph
      .get_nodes()
      .map(|node| tree.names[node.key()].to_owned())
      .sorted()
      .collect_vec();
    let internal_names: BTreeSet<String> = tree
      .graph
      .get_internal_nodes()
      .map(|node| tree.names[node.key()].to_owned())
      .collect();
    let confidence_names = table_read_file::<NamedRow>(
      OutputSelection::ConfidenceTsv,
      &output.path().join("timetree.confidence.tsv"),
    )?
    .into_iter()
    .map(|row| row.name)
    .sorted()
    .collect_vec();
    let fasta_names: BTreeSet<String> =
      fasta_read_file(output.path().join("timetree.sequences.fasta"), &Alphabet::default())?
        .into_iter()
        .map(|record| record.seq_name)
        .collect();

    assert_eq!(BTreeSet::new(), tree.unnamed_nodes);
    assert_eq!(node_names, confidence_names);
    assert_eq!(internal_names, fasta_names);
    Ok(())
  }

  #[derive(Deserialize)]
  struct NamedRow {
    name: String,
  }
}
