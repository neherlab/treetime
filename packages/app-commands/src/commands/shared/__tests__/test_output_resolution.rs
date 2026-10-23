#[cfg(test)]
mod tests {
  use crate::commands::clock::args::TreetimeClockArgsRaw;
  use crate::commands::shared::output_args::{NewickAnnotationsArg, OutputCoreArgs, ResolveOutputs};
  use app_output::command_outputs::CommandKind;
  use app_output::output_conditions::OutputSettings;
  use app_output::output_plan::TreeWriteKind;
  use app_output::output_plan::{OutputSelection, PlannedFile, Requested, ResolvedOutputs};
  use clap::Parser;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::num::NonZeroUsize;
  use std::path::{Path, PathBuf};
  use tempfile::TempDir;
  use treetime_io::nwk::NwkStyle;
  use treetime_utils::assert_error;

  #[test]
  fn test_resolve_output_all_default_selection() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();

    assert!(helpers::tree_paths(&resolved).contains_key(&TreeWriteKind::Nwk(NwkStyle::Plain)));
    assert!(helpers::tree_paths(&resolved).contains_key(&TreeWriteKind::Nexus));
    assert!(
      !helpers::tree_paths(&resolved).contains_key(&TreeWriteKind::GraphJson),
      "GraphJson is not a default output"
    );
    assert!(
      !helpers::tree_paths(&resolved).contains_key(&TreeWriteKind::Dot),
      "Dot is not a default output"
    );
  }

  #[test]
  fn test_resolve_output_all_with_selection_filter() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[OutputSelection::Nwk],
      &[],
    )
    .unwrap();

    assert_eq!(helpers::tree_paths(&resolved).len(), 1);
    assert!(helpers::tree_paths(&resolved).contains_key(&TreeWriteKind::Nwk(NwkStyle::Plain)));
  }

  #[test]
  fn test_resolve_per_file_only_no_output_all() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("/tmp/test.nwk")),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();

    assert_eq!(
      helpers::tree_paths(&resolved),
      btreemap! { TreeWriteKind::Nwk(NwkStyle::Plain) => PathBuf::from("/tmp/test.nwk") }
    );
  }

  #[test]
  fn test_resolve_no_outputs_errors() {
    let args = OutputCoreArgs::default();
    let result = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    );
    assert_error!(
      result,
      "No outputs selected: add -O <DIR> to write the default outputs into a folder, or a per-file --output-* flag"
    );
  }

  #[test]
  fn test_resolve_rejects_duplicate_output_destinations() {
    let path = PathBuf::from("same-output");
    let args = OutputCoreArgs {
      output_tree_nwk: Some(path.clone()),
      ..OutputCoreArgs::default()
    };

    let result = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[(OutputSelection::Model, Some(path.as_path()))],
    );

    assert_error!(
      result,
      "Output destination 'same-output' is selected more than once (--output-tree-nwk and --output-model)"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::dash_then_dev_stdout( ("-",           "/dev/stdout"), "Output destination '/dev/stdout' is selected more than once (--output-tree-nwk and --output-model)")]
  #[case::dev_stdout_then_dash( ("/dev/stdout", "-"),           "Output destination '-' is selected more than once (--output-tree-nwk and --output-model)")]
  #[trace]
  fn test_resolve_rejects_two_outputs_to_stdout_whatever_the_spelling(
    #[case] (tree_path, gtr_path): (&str, &str),
    #[case] expected: &str,
  ) {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from(tree_path)),
      ..OutputCoreArgs::default()
    };

    let result = helpers::resolve(&args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[(OutputSelection::Model, Some(Path::new(gtr_path)))],
    );

    assert_error!(result, expected);
  }

  #[test]
  fn test_resolve_rejects_two_tree_outputs_with_one_destination() {
    let path = PathBuf::from("same-output");
    let args = OutputCoreArgs {
      output_tree_nwk: Some(path.clone()),
      output_tree_nexus: Some(path),
      ..OutputCoreArgs::default()
    };

    let result = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    );

    assert_error!(
      result,
      "Output destination 'same-output' is selected more than once (--output-tree-nwk and --output-tree-nexus)"
    );
  }

  #[test]
  fn test_resolve_non_tree_explicit_path() {
    let dir = TempDir::new().unwrap();
    let gtr_path = dir.path().join("my.model.json");
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[(OutputSelection::Model, Some(gtr_path.as_path()))],
    )
    .unwrap();

    let expected = PlannedFile {
      path: gtr_path,
      requested: Requested::Named,
      condition: None,
    };
    assert_eq!(expected, resolved.non_tree_outputs[&OutputSelection::Model]);
  }

  #[test]
  fn test_resolve_non_tree_default_path_uses_command_stem() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[(OutputSelection::Model, None)],
    )
    .unwrap();

    let expected = PlannedFile {
      path: dir.path().join("ancestral.model.json"),
      requested: Requested::All,
      condition: None,
    };
    assert_eq!(expected, resolved.non_tree_outputs[&OutputSelection::Model]);
  }

  #[test]
  fn test_resolve_tree_default_paths_use_command_stem() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Clock,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[OutputSelection::Nwk],
      &[],
    )
    .unwrap();
    assert_eq!(
      &helpers::tree_paths(&resolved)[&TreeWriteKind::Nwk(NwkStyle::Plain)],
      &dir.path().join("clock.nwk")
    );
  }

  #[test]
  fn test_resolve_timetree_defaults_include_auspice() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Timetree,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();
    assert!(helpers::tree_paths(&resolved).contains_key(&TreeWriteKind::Auspice));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral(CommandKind::Ancestral)]
  #[case::timetree (CommandKind::Timetree)]
  #[case::optimize (CommandKind::Optimize)]
  #[case::mugration(CommandKind::Mugration)]
  #[case::clock    (CommandKind::Clock)]
  #[case::prune    (CommandKind::Prune)]
  #[trace]
  fn test_resolve_all_expands_to_complete_tree_set(#[case] command: CommandKind) {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(&args, command, NonZeroUsize::MIN, OutputSettings::default(), &[OutputSelection::All], &[]).unwrap();
    let actual = helpers::tree_paths(&resolved).keys().copied().collect();
    let expected = btreeset! {
      TreeWriteKind::Nwk(NwkStyle::Plain),
      TreeWriteKind::Nexus,
      TreeWriteKind::Auspice,
      TreeWriteKind::MatPb,
      TreeWriteKind::MatJson,
      TreeWriteKind::GraphJson,
      TreeWriteKind::Dot,
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_resolve_explicit_mat_available_on_clock() {
    let args = OutputCoreArgs {
      output_tree_mat_pb: Some(PathBuf::from("/tmp/out.pb")),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Clock,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();
    assert_eq!(
      &PathBuf::from("/tmp/out.pb"),
      &helpers::tree_paths(&resolved)[&TreeWriteKind::MatPb]
    );
  }

  #[test]
  fn test_resolve_explicit_auspice_available_on_clock() {
    let args = OutputCoreArgs {
      output_tree_auspice: Some(PathBuf::from("/tmp/out.auspice.json")),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Clock,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();
    assert_eq!(
      &PathBuf::from("/tmp/out.auspice.json"),
      &helpers::tree_paths(&resolved)[&TreeWriteKind::Auspice]
    );
  }

  #[test]
  fn test_s2_row1_output_all_default_annotations() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[OutputSelection::Nwk, OutputSelection::Nexus],
      &[],
    )
    .unwrap();
    let expected = btreemap! {
      TreeWriteKind::Nwk(NwkStyle::Plain) => dir.path().join("ancestral.nwk"),
      TreeWriteKind::Nexus => dir.path().join("ancestral.nexus"),
    };
    assert_eq!(helpers::tree_paths(&resolved), expected);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast(NewickAnnotationsArg::Beast, NwkStyle::Beast)]
  #[case::nhx(  NewickAnnotationsArg::Nhx,   NwkStyle::Nhx)]
  #[trace]
  fn test_s2_row2_output_all_newick_annotations_keep_the_file_names(
    #[case] newick_annotations: NewickAnnotationsArg,
    #[case] style: NwkStyle,
  ) {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      newick_annotations,
      ..Default::default()
    };
    let resolved = helpers::resolve(&args,
        CommandKind::Ancestral,
        NonZeroUsize::MIN,
        OutputSettings::default(),
        &[OutputSelection::Nwk, OutputSelection::Nexus],
        &[],
      )
      .unwrap();
    let expected = btreemap! {
      TreeWriteKind::Nwk(style) => dir.path().join("ancestral.nwk"),
      TreeWriteKind::Nexus => dir.path().join("ancestral.nexus"),
    };
    assert_eq!(helpers::tree_paths(&resolved), expected);
  }

  #[test]
  fn test_s2_row3_per_file_default_annotations() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();
    assert_eq!(
      helpers::tree_paths(&resolved),
      btreemap! { TreeWriteKind::Nwk(NwkStyle::Plain) => PathBuf::from("my.nwk") }
    );
  }

  #[test]
  fn test_s2_row4_per_file_with_annotations_used_as_is() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      newick_annotations: NewickAnnotationsArg::Beast,
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();
    assert_eq!(
      helpers::tree_paths(&resolved),
      btreemap! { TreeWriteKind::Nwk(NwkStyle::Beast) => PathBuf::from("my.nwk") }
    );
  }

  #[test]
  fn test_s2_row5_per_file_overrides_output_all_nwk_only() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[OutputSelection::Nwk, OutputSelection::Nexus],
      &[],
    )
    .unwrap();
    let expected = btreemap! {
      TreeWriteKind::Nwk(NwkStyle::Plain) => PathBuf::from("my.nwk"),
      TreeWriteKind::Nexus => dir.path().join("ancestral.nexus"),
    };
    assert_eq!(helpers::tree_paths(&resolved), expected);
  }

  #[test]
  fn test_s2_row6_per_file_supplements_selection() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[OutputSelection::Nexus],
      &[],
    )
    .unwrap();
    let expected = btreemap! {
      TreeWriteKind::Nwk(NwkStyle::Plain) => PathBuf::from("my.nwk"),
      TreeWriteKind::Nexus => dir.path().join("ancestral.nexus"),
    };
    assert_eq!(helpers::tree_paths(&resolved), expected);
  }

  #[test]
  fn test_s2_row7_newick_annotations_leave_the_nexus_comments_unchanged() {
    let args = OutputCoreArgs {
      output_tree_nwk: Some(PathBuf::from("my.nwk")),
      output_tree_nexus: Some(PathBuf::from("my.nexus")),
      newick_annotations: NewickAnnotationsArg::Nhx,
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Ancestral,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[],
      &[],
    )
    .unwrap();
    let expected: BTreeMap<TreeWriteKind, PathBuf> = btreemap! {
      TreeWriteKind::Nwk(NwkStyle::Nhx) => PathBuf::from("my.nwk"),
      TreeWriteKind::Nexus => PathBuf::from("my.nexus"),
    };
    assert_eq!(helpers::tree_paths(&resolved), expected);
  }

  #[test]
  fn test_resolve_timetree_default_with_a_prior_includes_the_coalescent_table_not_json() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Timetree,
      NonZeroUsize::MIN,
      OutputSettings {
        coalescent: true,
        ..OutputSettings::default()
      },
      &[],
      &[],
    )
    .unwrap();

    assert_eq!(
      Some(dir.path().join("timetree.coalescent.tsv").as_path()),
      resolved.path(OutputSelection::Coalescent)
    );
    assert!(
      !resolved.non_tree_outputs.contains_key(&OutputSelection::CoalescentJson),
      "coalescent JSON is opt-in, not a default output"
    );
  }

  #[test]
  fn test_resolve_timetree_coalescent_json_opt_in_via_selection() {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };
    let resolved = helpers::resolve(
      &args,
      CommandKind::Timetree,
      NonZeroUsize::MIN,
      OutputSettings::default(),
      &[OutputSelection::CoalescentJson],
      &[],
    )
    .unwrap();

    let expected = btreemap! {
      OutputSelection::CoalescentJson => PlannedFile {
        path: dir.path().join("timetree.coalescent.json"),
        requested: Requested::All,
        condition: None,
      },
    };
    assert_eq!(expected, resolved.non_tree_outputs);
  }

  #[test]
  fn test_resolve_timetree_offers_both_coalescent_outputs() {
    let coalescent = btreeset! {
      OutputSelection::Coalescent,
      OutputSelection::CoalescentJson,
    };
    assert!(coalescent.is_subset(&CommandKind::Timetree.all_selectable(OutputSettings::default())));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ancestral(CommandKind::Ancestral)]
  #[case::optimize (CommandKind::Optimize)]
  #[case::mugration(CommandKind::Mugration)]
  #[case::clock    (CommandKind::Clock)]
  #[case::prune    (CommandKind::Prune)]
  #[trace]
  fn test_resolve_coalescent_outputs_are_timetree_only(#[case] command: CommandKind) {
    let coalescent = btreeset! {
      OutputSelection::Coalescent,
      OutputSelection::CoalescentJson,
    };
    assert!(
      command.all_selectable(OutputSettings::default()).is_disjoint(&coalescent),
      "{command:?} must not offer coalescent outputs"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_set(   vec![],                          vec![("clock.rtt.svg", OutputSelection::RttChart)])]
  #[case::nwk_only(      vec![OutputSelection::Nwk],      vec![])]
  #[case::chart_selected(vec![OutputSelection::RttChart], vec![("clock.rtt.svg", OutputSelection::RttChart)])]
  #[trace]
  fn test_resolve_clock_charts_follow_the_selection(
    #[case] selection: Vec<OutputSelection>,
    #[case] expected: Vec<(&str, OutputSelection)>,
  ) {
    let dir = TempDir::new().unwrap();
    let args = OutputCoreArgs {
      output_all: Some(dir.path().to_path_buf()),
      ..Default::default()
    };

    let resolved = helpers::resolve(&args, CommandKind::Clock, NonZeroUsize::MIN, OutputSettings::default(), &selection, &[]).unwrap();

    let expected: BTreeMap<OutputSelection, PlannedFile> = expected
      .into_iter()
      .map(|(name, selection)| {
        let file = PlannedFile {
          path: dir.path().join(name),
          requested: Requested::All,
          condition: None,
        };
        (selection, file)
      })
      .collect();
    let actual: BTreeMap<OutputSelection, PlannedFile> = resolved
      .non_tree_outputs
      .into_iter()
      .filter(|(selection, _)| *selection == OutputSelection::RttChart)
      .collect();
    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_set(      &["--output-all={dir}"],                                             &[("clock.rtt.svg", Requested::All)])]
  #[case::pdf_path(         &["--output-rtt-chart={dir}/chart.pdf"],                             &[("chart.pdf",     Requested::Named)])]
  #[case::pdf_path_over_all(&["--output-all={dir}", "--output-rtt-chart={dir}/chart.pdf"],       &[("chart.pdf",     Requested::Named)])]
  #[trace]
  fn test_resolve_rtt_chart_from_command_line(#[case] flags: &[&str], #[case] expected: &[(&str, Requested)]) {
    let dir = TempDir::new().unwrap();
    let dir_str = dir.path().to_str().unwrap();
    let args = TreetimeClockArgsRaw::try_parse_from(
      ["treetime", "--metadata=/dev/null"]
        .into_iter()
        .map(str::to_owned)
        .chain(flags.iter().map(|flag| flag.replace("{dir}", dir_str))),
    )
    .unwrap();

    let resolved = args.resolve_outputs().unwrap();

    let expected: BTreeMap<OutputSelection, PlannedFile> = expected
      .iter()
      .map(|&(name, requested)| {
        let file = PlannedFile {
          path: dir.path().join(name),
          requested,
        };
        (OutputSelection::RttChart, file)
      })
      .collect();
    let actual: BTreeMap<OutputSelection, PlannedFile> = resolved
      .non_tree_outputs
      .into_iter()
      .filter(|(selection, _)| *selection == OutputSelection::RttChart)
      .collect();
    assert_eq!(expected, actual);
  }

  mod helpers {
    use super::*;
    use crate::commands::shared::output_args::{OutputRequest, TreeOutputPaths};
    use eyre::Report;

    pub(super) fn resolve(
      args: &OutputCoreArgs,
      command: CommandKind,
      alignment_count: NonZeroUsize,
      settings: OutputSettings,
      selection: &[OutputSelection],
      file_paths: &[(OutputSelection, Option<&Path>)],
    ) -> Result<ResolvedOutputs, Report> {
      let args = args.clone();
      OutputRequest {
        command,
        alignment_count: alignment_count.get(),
        settings,
        output_all: args.output_all,
        newick_annotations: args.newick_annotations,
        selection: selection.to_vec(),
        tree_paths: TreeOutputPaths {
          nwk: args.output_tree_nwk,
          nexus: args.output_tree_nexus,
          auspice: args.output_tree_auspice,
          mat_pb: args.output_tree_mat_pb,
          mat_json: args.output_tree_mat_json,
          graph_json: args.output_tree_graph_json,
          dot: args.output_tree_dot,
        },
        file_paths: file_paths
          .iter()
          .map(|(selection, path)| (*selection, path.map(Path::to_path_buf)))
          .collect(),
      }
      .resolve()
    }

    pub(super) fn tree_paths(resolved: &ResolvedOutputs) -> BTreeMap<TreeWriteKind, PathBuf> {
      let mut paths = resolved.alignment_tree_outputs(1);
      paths.extend(resolved.tree_outputs.clone());
      paths
    }
  }
}
