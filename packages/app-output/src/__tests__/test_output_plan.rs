#[cfg(test)]
mod tests {
  use crate::output_plan::{
    CommandKind, OutputSelection, PlannedFile, Requested, TreeWriteKind, output_unavailable, plan,
  };
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::{Path, PathBuf};
  use std::str::FromStr;
  use strum::IntoEnumIterator;
  use treetime::progress::LogLevel;
  use treetime_io::nwk::NwkStyle;
  use treetime_utils::io::json::to_json_value;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_output_plan_selection_tag_matches_serde_name() {
    let serde_names = OutputSelection::iter()
      .map(|selection| to_json_value(&selection).unwrap())
      .collect_vec();
    let tag_names = OutputSelection::iter()
      .map(|selection| serde_json::Value::String(selection.as_ref().to_owned()))
      .collect_vec();
    assert_eq!(serde_names, tag_names);
  }

  #[test]
  fn test_output_plan_selection_from_str_roundtrip() {
    let expected = OutputSelection::iter().collect_vec();
    let actual = OutputSelection::iter()
      .map(|selection| OutputSelection::from_str(selection.as_ref()).unwrap())
      .collect_vec();
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_output_plan_selection_from_str_rejects_unknown_name() {
    assert_eq!(
      Err(strum::ParseError::VariantNotFound),
      OutputSelection::from_str("mat_pb")
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::default_set(  vec![],                                                    Requested::All)]
  #[case::listed_value( vec![OutputSelection::Model],                                Requested::All)]
  #[case::all_value(    vec![OutputSelection::All],                                Requested::All)]
  #[trace]
  fn test_output_plan_output_all_entries_are_requested_by_output_all(
    #[case] selection: Vec<OutputSelection>,
    #[case] expected: Requested,
  ) -> Result<(), Report> {
    let request = helpers::request(selection, btreemap! {});

    let resolved = plan(&request)?;

    assert_eq!(expected, resolved.non_tree_outputs[&OutputSelection::Model].requested);
    Ok(())
  }

  #[test]
  fn test_output_plan_per_file_flag_is_named_even_with_output_all() -> Result<(), Report> {
    let request = helpers::request(
      vec![],
      btreemap! { OutputSelection::Model => PathBuf::from("my.model.json") },
    );

    let resolved = plan(&request)?;

    let expected = PlannedFile {
      path: PathBuf::from("my.model.json"),
      requested: Requested::Named,
    };
    assert_eq!(expected, resolved.non_tree_outputs[&OutputSelection::Model]);
    Ok(())
  }

  #[test]
  fn test_output_plan_unavailable_named_output_fails_with_flag_and_reason() {
    let file = PlannedFile {
      path: PathBuf::from("my.model.json"),
      requested: Requested::Named,
    };
    let (log, messages) = helpers::recording_log();

    assert_error!(
      output_unavailable(OutputSelection::Model, &file, "no GTR model was fitted", &log),
      "--output-model was requested, but no GTR model was fitted"
    );
    drop(log);
    assert_eq!(Vec::<(LogLevel, String)>::new(), messages.iter().collect_vec());
  }

  #[test]
  fn test_output_plan_unavailable_output_from_output_all_is_skipped_with_debug_message() -> Result<(), Report> {
    let file = PlannedFile {
      path: PathBuf::from("out/ancestral.model.json"),
      requested: Requested::All,
    };
    let (log, messages) = helpers::recording_log();

    output_unavailable(OutputSelection::Model, &file, "no GTR model was fitted", &log)?;

    drop(log);
    let expected = vec![(
      LogLevel::Debug,
      o!("Not writing 'out/ancestral.model.json': no GTR model was fitted"),
    )];
    assert_eq!(expected, messages.iter().collect_vec());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::two_parts(       ("out/aa.fasta",     "out/aa.fasta"), "Output destination 'out/aa.fasta' is selected more than once (part S and part M)")]
  #[case::stdout_spellings(("-",                "/dev/stdout"),  "Output destination '/dev/stdout' is selected more than once (part S and part M)")]
  #[case::planned_output(  ("out/gtr.json",     "out/M.fasta"),  "Output destination 'out/gtr.json' is selected more than once (--output-model and part S)")]
  #[trace]
  fn test_output_plan_expansion_rejects_shared_destination(
    #[case] (path_s, path_m): (&str, &str),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let resolved = plan(&helpers::request(
      vec![OutputSelection::Model],
      btreemap! {
        OutputSelection::Model => PathBuf::from("out/gtr.json"),
        OutputSelection::Translations => PathBuf::from("out/aa.fasta"),
      },
    ))?;

    assert_error!(
      resolved.ensure_unique_with_expansion(
        OutputSelection::Translations,
        [(o!("part S"), Path::new(path_s)), (o!("part M"), Path::new(path_m))],
      ),
      expected
    );
    Ok(())
  }

  #[test]
  fn test_output_plan_expansion_replaces_the_expanded_output() -> Result<(), Report> {
    let resolved = plan(&helpers::request(
      vec![OutputSelection::Model],
      btreemap! { OutputSelection::Translations => PathBuf::from("out/aa.fasta") },
    ))?;

    resolved.ensure_unique_with_expansion(
      OutputSelection::Translations,
      [(o!("part S"), Path::new("out/aa.fasta"))],
    )
  }

  #[test]
  fn test_output_plan_homoplasy_default_outputs() -> Result<(), Report> {
    let resolved = plan(&helpers::homoplasy_request(btreemap! {}))?;

    let expected = btreemap! {
      OutputSelection::Nwk => vec![PathBuf::from("out/homoplasy.nwk")],
      OutputSelection::Nexus => vec![PathBuf::from("out/homoplasy.nexus")],
      OutputSelection::HomoplasyStats => vec![PathBuf::from("out/homoplasy.stats.json")],
      OutputSelection::HomoplasyReport => vec![PathBuf::from("out/homoplasy.report.txt")],
    };
    assert_eq!(expected, resolved.paths_by_selection());
    Ok(())
  }

  #[test]
  fn test_output_plan_homoplasy_per_file_flags_override_output_all() -> Result<(), Report> {
    let resolved = plan(&helpers::homoplasy_request(btreemap! {
      OutputSelection::HomoplasyStats => PathBuf::from("stats.json"),
      OutputSelection::HomoplasyReport => PathBuf::from("-"),
    }))?;

    assert_eq!(
      (Some(Path::new("stats.json")), Some(Path::new("-"))),
      (
        resolved.path(OutputSelection::HomoplasyStats),
        resolved.path(OutputSelection::HomoplasyReport)
      )
    );
    Ok(())
  }

  #[test]
  fn test_output_plan_one_partition_keeps_default_names() -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      1,
      Some("out"),
      btreemap! {},
      btreemap! {},
    ))?;

    let expected = btreemap! {
      OutputSelection::Nwk => vec![PathBuf::from("out/ancestral.nwk")],
      OutputSelection::Nexus => vec![PathBuf::from("out/ancestral.nexus")],
      OutputSelection::Auspice => vec![PathBuf::from("out/ancestral.auspice.json")],
      OutputSelection::MatPb => vec![PathBuf::from("out/ancestral.mat.pb")],
      OutputSelection::NodeData => vec![PathBuf::from("out/ancestral.node-data.json")],
      OutputSelection::Model => vec![PathBuf::from("out/ancestral.model.json")],
      OutputSelection::Sequences => vec![PathBuf::from("out/ancestral.sequences.fasta")],
    };
    assert_eq!(expected, resolved.paths_by_selection());
    Ok(())
  }

  #[test]
  fn test_output_plan_two_partitions_number_default_names_after_the_stem() -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      2,
      Some("out"),
      btreemap! {},
      btreemap! {},
    ))?;

    let expected = btreemap! {
      OutputSelection::Nwk => vec![PathBuf::from("out/ancestral.nwk")],
      OutputSelection::Nexus => vec![PathBuf::from("out/ancestral.1.nexus"), PathBuf::from("out/ancestral.2.nexus")],
      OutputSelection::Auspice => vec![
        PathBuf::from("out/ancestral.1.auspice.json"),
        PathBuf::from("out/ancestral.2.auspice.json"),
      ],
      OutputSelection::MatPb => vec![PathBuf::from("out/ancestral.1.mat.pb"), PathBuf::from("out/ancestral.2.mat.pb")],
      OutputSelection::NodeData => vec![
        PathBuf::from("out/ancestral.1.node-data.json"),
        PathBuf::from("out/ancestral.2.node-data.json"),
      ],
      OutputSelection::Model => vec![PathBuf::from("out/ancestral.model.json")],
      OutputSelection::Sequences => vec![
        PathBuf::from("out/ancestral.1.sequences.fasta"),
        PathBuf::from("out/ancestral.2.sequences.fasta"),
      ],
    };
    assert_eq!(expected, resolved.paths_by_selection());
    Ok(())
  }

  #[test]
  fn test_output_plan_partition_lookups_return_the_path_of_each_partition() -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      2,
      Some("out"),
      btreemap! {},
      btreemap! {},
    ))?;

    let actual = (
      [1, 2, 3].map(|partition| resolved.partition_path(OutputSelection::Sequences, partition)),
      [0, 2].map(|partition| resolved.tree_partition_path(TreeWriteKind::Nexus, partition)),
      resolved.tree_partition_path(TreeWriteKind::Nwk(NwkStyle::Plain), 1),
      resolved.partition_path(OutputSelection::Model, 1),
    );

    let expected = (
      [
        Some(Path::new("out/ancestral.1.sequences.fasta")),
        Some(Path::new("out/ancestral.2.sequences.fasta")),
        None,
      ],
      [None, Some(Path::new("out/ancestral.2.nexus"))],
      None,
      None,
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_partition( 1)]
  #[case::two_partitions(2)]
  #[trace]
  fn test_output_plan_per_partition_outputs_have_no_run_path(#[case] partition_count: usize) -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      partition_count,
      Some("out"),
      btreemap! {},
      btreemap! {},
    ))?;

    let actual = (
      resolved.path(OutputSelection::Sequences),
      resolved.path(OutputSelection::NodeData),
      resolved.path(OutputSelection::Model),
      resolved.tree_outputs.keys().copied().collect_vec(),
      resolved.partition_tree_outputs(partition_count).into_keys().collect_vec(),
    );

    let expected = (
      None,
      None,
      Some(Path::new("out/ancestral.model.json")),
      vec![TreeWriteKind::Nwk(NwkStyle::Plain)],
      vec![TreeWriteKind::Nexus, TreeWriteKind::Auspice, TreeWriteKind::MatPb],
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_output_plan_trees_of_commands_without_alignments_are_written_once() -> Result<(), Report> {
    let mut request = helpers::partitioned_request(1, Some("out"), btreemap! {}, btreemap! {});
    request.command = CommandKind::Mugration;
    request.selection = vec![
      OutputSelection::Nwk,
      OutputSelection::Nexus,
      OutputSelection::Auspice,
      OutputSelection::MatPb,
    ];
    request.newick_annotations = NwkStyle::Beast;

    let resolved = plan(&request)?;

    let expected = (
      vec![
        TreeWriteKind::Nwk(NwkStyle::Beast),
        TreeWriteKind::Nexus,
        TreeWriteKind::Auspice,
        TreeWriteKind::MatPb,
      ],
      true,
    );
    let actual = (
      resolved.tree_outputs.keys().copied().collect_vec(),
      resolved.tree_outputs_by_partition.is_empty(),
    );
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_output_plan_tree_of_a_command_without_alignments_rejects_the_placeholder() {
    let mut request = helpers::partitioned_request(
      1,
      None,
      btreemap! { OutputSelection::Auspice => PathBuf::from("tree.{alignment}.json") },
      btreemap! {},
    );
    request.command = CommandKind::Clock;

    assert_error!(
      plan(&request),
      "--output-tree-auspice is written once per run, so its file name cannot contain '{alignment}'"
    );
  }

  #[test]
  fn test_output_plan_explicit_paths_expand_the_partition_placeholder() -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      2,
      None,
      btreemap! {
        OutputSelection::Nwk => PathBuf::from("trees/{alignment}.nwk"),
        OutputSelection::Nexus => PathBuf::from("trees/{alignment}.nexus"),
        OutputSelection::Auspice => PathBuf::from("auspice/segment-{alignment}.json"),
        OutputSelection::GraphJson => PathBuf::from("graph.json"),
        OutputSelection::Dot => PathBuf::from("graph.{alignment}.dot"),
      },
      btreemap! {
        OutputSelection::Sequences => PathBuf::from("seq.{alignment}.fasta"),
        OutputSelection::Model => PathBuf::from("gtr.json"),
      },
    ))?;

    let expected = btreemap! {
      OutputSelection::Nwk => vec![PathBuf::from("trees/1.nwk"), PathBuf::from("trees/2.nwk")],
      OutputSelection::Nexus => vec![PathBuf::from("trees/1.nexus"), PathBuf::from("trees/2.nexus")],
      OutputSelection::Auspice => vec![
        PathBuf::from("auspice/segment-1.json"),
        PathBuf::from("auspice/segment-2.json"),
      ],
      OutputSelection::GraphJson => vec![PathBuf::from("graph.json")],
      OutputSelection::Dot => vec![PathBuf::from("graph.1.dot"), PathBuf::from("graph.2.dot")],
      OutputSelection::Model => vec![PathBuf::from("gtr.json")],
      OutputSelection::Sequences => vec![PathBuf::from("seq.1.fasta"), PathBuf::from("seq.2.fasta")],
    };
    assert_eq!(expected, resolved.paths_by_selection());
    Ok(())
  }

  #[test]
  fn test_output_plan_one_partition_expands_the_placeholder_to_one() -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      1,
      None,
      btreemap! {},
      btreemap! { OutputSelection::Sequences => PathBuf::from("seq.{alignment}.fasta") },
    ))?;

    assert_eq!(
      Some(Path::new("seq.1.fasta")),
      resolved.partition_path(OutputSelection::Sequences, 1)
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::fasta(   OutputSelection::Sequences, "seq.fasta",    "--output-sequences")]
  #[case::node_data(OutputSelection::NodeData,        "node.json",    "--output-node-data")]
  #[case::stdout(  OutputSelection::Sequences, "-",            "--output-sequences")]
  #[case::folder(  OutputSelection::Sequences, "{alignment}/seq.fasta", "--output-sequences")]
  #[trace]
  fn test_output_plan_two_partitions_need_the_placeholder_in_explicit_paths(
    #[case] selection: OutputSelection,
    #[case] path: &str,
    #[case] flag: &str,
  ) {
    assert_error!(
      plan(&helpers::partitioned_request(2, None, btreemap! {}, btreemap! { selection => PathBuf::from(path) })),
      format!(
        "{flag} needs the placeholder '{{alignment}}' in its file name when the run has several \
         alignment files, otherwise the output of each alignment file overwrites the same file"
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nexus(        OutputSelection::Nexus, "tree.nexus", NwkStyle::Plain, "--output-tree-nexus")]
  #[case::newick_beast( OutputSelection::Nwk,   "tree.nwk",   NwkStyle::Beast, "--output-tree-nwk")]
  #[case::newick_nhx(   OutputSelection::Nwk,   "tree.nwk",   NwkStyle::Nhx,   "--output-tree-nwk")]
  #[trace]
  fn test_output_plan_two_partitions_need_the_placeholder_in_tree_paths_with_comments(
    #[case] selection: OutputSelection,
    #[case] path: &str,
    #[case] newick_annotations: NwkStyle,
    #[case] flag: &str,
  ) {
    let mut request = helpers::partitioned_request(2, None, btreemap! { selection => PathBuf::from(path) }, btreemap! {});
    request.newick_annotations = newick_annotations;
    assert_error!(
      plan(&request),
      format!(
        "{flag} needs the placeholder '{{alignment}}' in its file name when the run has several \
         alignment files, otherwise the output of each alignment file overwrites the same file"
      )
    );
  }

  #[test]
  fn test_output_plan_two_partitions_write_plain_newick_once() -> Result<(), Report> {
    let resolved = plan(&helpers::partitioned_request(
      2,
      None,
      btreemap! { OutputSelection::Nwk => PathBuf::from("tree.nwk") },
      btreemap! {},
    ))?;

    let expected = btreemap! { TreeWriteKind::Nwk(NwkStyle::Plain) => PathBuf::from("tree.nwk") };
    assert_eq!(expected, resolved.tree_outputs);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast(NwkStyle::Beast)]
  #[case::nhx(  NwkStyle::Nhx)]
  #[trace]
  fn test_output_plan_two_partitions_write_annotated_newick_once_per_alignment(
    #[case] newick_annotations: NwkStyle,
  ) -> Result<(), Report> {
    let mut request = helpers::partitioned_request(2, Some("out"), btreemap! {}, btreemap! {});
    request.selection = vec![OutputSelection::Nwk];
    request.newick_annotations = newick_annotations;

    let resolved = plan(&request)?;

    let expected = btreemap! {
      TreeWriteKind::Nwk(newick_annotations) => vec![PathBuf::from("out/ancestral.1.nwk"), PathBuf::from("out/ancestral.2.nwk")],
    };
    assert_eq!(expected, resolved.tree_outputs_by_partition);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_partition( 1)]
  #[case::two_partitions(2)]
  #[trace]
  fn test_output_plan_single_file_output_rejects_the_placeholder(#[case] partition_count: usize) {
    assert_error!(
      plan(&helpers::partitioned_request(
        partition_count,
        Some("out"),
        btreemap! {},
        btreemap! { OutputSelection::Model => PathBuf::from("gtr.{alignment}.json") },
      )),
      "--output-model is written once per run, so its file name cannot contain '{alignment}'"
    );
  }

  #[test]
  fn test_output_plan_partition_path_colliding_with_another_output_fails() {
    assert_error!(
      plan(&helpers::partitioned_request(
        2,
        None,
        btreemap! {},
        btreemap! {
          OutputSelection::Sequences => PathBuf::from("out.{alignment}.fasta"),
          OutputSelection::Model => PathBuf::from("out.2.fasta"),
        },
      )),
      "Output destination 'out.2.fasta' is selected more than once (--output-model and --output-sequences for alignment 2)"
    );
  }

  #[test]
  fn test_output_plan_timetree_node_data_stays_one_file_with_two_partitions() -> Result<(), Report> {
    let mut request = helpers::partitioned_request(2, Some("out"), btreemap! {}, btreemap! {});
    request.command = CommandKind::Timetree;
    request.selection = vec![OutputSelection::NodeData];

    let resolved = plan(&request)?;

    let expected = btreemap! {
      OutputSelection::NodeData => vec![PathBuf::from("out/timetree.node-data.json")],
    };
    assert_eq!(expected, resolved.paths_by_selection());
    Ok(())
  }

  mod helpers {
    use crate::output_plan::{CommandKind, OutputPlanRequest, OutputSelection};
    use std::collections::BTreeMap;
    use std::num::NonZeroUsize;
    use std::path::PathBuf;
    use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
    use treetime::progress::{LogLevel, LogSink};
    use treetime_io::nwk::NwkStyle;

    pub(super) fn request(
      selection: Vec<OutputSelection>,
      non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
    ) -> OutputPlanRequest {
      OutputPlanRequest {
        command: CommandKind::Ancestral,
        partition_count: NonZeroUsize::MIN,
        output_all: Some(PathBuf::from("out")),
        newick_annotations: NwkStyle::Plain,
        selection,
        tree_overrides: BTreeMap::new(),
        non_tree_overrides,
      }
    }

    pub(super) fn partitioned_request(
      partition_count: usize,
      output_all: Option<&str>,
      tree_overrides: BTreeMap<OutputSelection, PathBuf>,
      non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
    ) -> OutputPlanRequest {
      OutputPlanRequest {
        command: CommandKind::Ancestral,
        partition_count: NonZeroUsize::new(partition_count).expect("tests plan at least one partition"),
        selection: if output_all.is_some() {
          vec![
            OutputSelection::Nwk,
            OutputSelection::Nexus,
            OutputSelection::Auspice,
            OutputSelection::MatPb,
            OutputSelection::NodeData,
            OutputSelection::Model,
            OutputSelection::Sequences,
          ]
        } else {
          vec![]
        },
        output_all: output_all.map(PathBuf::from),
        newick_annotations: NwkStyle::Plain,
        tree_overrides,
        non_tree_overrides,
      }
    }

    pub(super) fn homoplasy_request(non_tree_overrides: BTreeMap<OutputSelection, PathBuf>) -> OutputPlanRequest {
      OutputPlanRequest {
        command: CommandKind::Homoplasy,
        partition_count: NonZeroUsize::MIN,
        output_all: Some(PathBuf::from("out")),
        newick_annotations: NwkStyle::Plain,
        selection: vec![],
        tree_overrides: BTreeMap::new(),
        non_tree_overrides,
      }
    }

    pub(super) fn recording_log() -> (RecordingLog, Receiver<(LogLevel, String)>) {
      let (sender, receiver) = sync_channel(8);
      (RecordingLog { sender }, receiver)
    }

    pub(super) struct RecordingLog {
      sender: SyncSender<(LogLevel, String)>,
    }

    impl LogSink for RecordingLog {
      fn log(&self, level: LogLevel, message: &str) {
        self
          .sender
          .send((level, message.to_owned()))
          .expect("the test keeps the receiver alive");
      }

      fn log_enabled(&self, _level: LogLevel) -> bool {
        true
      }
    }
  }
}
