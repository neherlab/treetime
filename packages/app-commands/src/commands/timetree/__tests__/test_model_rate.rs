#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::project_root;
  use crate::commands::shared::alignment::AlignmentArgs;
  use crate::commands::shared::output_args::{OutputCoreArgs, TimetreeOutputSelection};
  use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
  use crate::commands::timetree::run::run_timetree_estimation;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::cancel::NoopCancel;
  use treetime::clock::clock_model::ClockModel;
  use treetime::gtr::get_gtr::{ModelEntry, ModelOutput};
  use treetime::progress::NoopProgress;
  use treetime_utils::io::json::json_read_file;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rustfmt::skip]
  #[rstest]
  #[case::dense( true )]
  #[case::sparse(false)]
  #[trace]
  fn test_timetree_one_alignment_with_inferred_model_has_unit_overall_rate(#[case] dense: bool) -> Result<(), Report> {
    let output = tempfile::tempdir()?;
    let dir = project_root().join("data/flu/h3n2/20");
    let args = TreetimeTimetreeArgs::try_from(TreetimeTimetreeArgsRaw {
      alignment: AlignmentArgs {
        alignment: vec![dir.join("aln.fasta.xz")],
      },
      tree: Some(dir.join("tree.nwk")),
      metadata: Some(dir.join("metadata.tsv")),
      dense: Some(dense),
      max_iter: 1,
      output: OutputCoreArgs {
        output_all: Some(output.path().to_path_buf()),
        ..OutputCoreArgs::default()
      },
      output_selection: vec![TimetreeOutputSelection::Gtr, TimetreeOutputSelection::ClockModel],
      ..TreetimeTimetreeArgsRaw::default()
    })?;

    run_timetree_estimation(&args, &NoopCancel, &NoopProgress, &NoopProgress)?;

    let models: ModelOutput = json_read_file(output.path().join("timetree.gtr.json"))?;
    let clock_model: ClockModel = json_read_file(output.path().join("timetree.clock-model.json"))?;
    let [ModelEntry::Sequence(model)] = models.models.as_slice() else {
      panic!("one alignment gives one sequence model");
    };
    pretty_assert_ulps_eq!(
      1.0,
      model.model.mu(),
      max_ulps = 0,
      "v0 sets the overall rate of an inferred model to 1: packages/legacy/treetime/treetime/treeanc.py#L1626-L1631"
    );
    assert_eq!(Some(clock_model.clock_rate()), model.clock_rate);
    Ok(())
  }
}
