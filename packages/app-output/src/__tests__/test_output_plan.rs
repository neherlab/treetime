#[cfg(test)]
mod tests {
  use crate::output_plan::{OutputSelection, PlannedFile, Requested, output_unavailable, plan};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::path::{Path, PathBuf};
  use std::str::FromStr;
  use strum::IntoEnumIterator;
  use treetime::progress::LogLevel;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_output_plan_selection_tag_matches_serde_name() {
    let serde_names = OutputSelection::iter()
      .map(|selection| serde_json::to_value(selection).unwrap())
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
  #[case::listed_value( vec![OutputSelection::Gtr],                                Requested::All)]
  #[case::all_value(    vec![OutputSelection::All],                                Requested::All)]
  #[trace]
  fn test_output_plan_output_all_entries_are_requested_by_output_all(
    #[case] selection: Vec<OutputSelection>,
    #[case] expected: Requested,
  ) -> Result<(), Report> {
    let request = helpers::request(selection, btreemap! {});

    let resolved = plan(&request)?;

    assert_eq!(expected, resolved.non_tree_outputs[&OutputSelection::Gtr].requested);
    Ok(())
  }

  #[test]
  fn test_output_plan_per_file_flag_is_named_even_with_output_all() -> Result<(), Report> {
    let request = helpers::request(
      vec![],
      btreemap! { OutputSelection::Gtr => PathBuf::from("my.gtr.json") },
    );

    let resolved = plan(&request)?;

    let expected = PlannedFile {
      path: PathBuf::from("my.gtr.json"),
      requested: Requested::Named,
    };
    assert_eq!(expected, resolved.non_tree_outputs[&OutputSelection::Gtr]);
    Ok(())
  }

  #[test]
  fn test_output_plan_unavailable_named_output_fails_with_flag_and_reason() {
    let file = PlannedFile {
      path: PathBuf::from("my.gtr.json"),
      requested: Requested::Named,
    };
    let (log, messages) = helpers::recording_log();

    assert_error!(
      output_unavailable(OutputSelection::Gtr, &file, "no GTR model was fitted", &log),
      "--output-gtr was requested, but no GTR model was fitted"
    );
    drop(log);
    assert_eq!(Vec::<(LogLevel, String)>::new(), messages.iter().collect_vec());
  }

  #[test]
  fn test_output_plan_unavailable_output_from_output_all_is_skipped_with_debug_message() -> Result<(), Report> {
    let file = PlannedFile {
      path: PathBuf::from("out/ancestral.gtr.json"),
      requested: Requested::All,
    };
    let (log, messages) = helpers::recording_log();

    output_unavailable(OutputSelection::Gtr, &file, "no GTR model was fitted", &log)?;

    drop(log);
    let expected = vec![(
      LogLevel::Debug,
      o!("Not writing 'out/ancestral.gtr.json': no GTR model was fitted"),
    )];
    assert_eq!(expected, messages.iter().collect_vec());
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::two_parts(       ("out/aa.fasta",     "out/aa.fasta"), "Output destination 'out/aa.fasta' is selected more than once (part S and part M)")]
  #[case::stdout_spellings(("-",                "/dev/stdout"),  "Output destination '/dev/stdout' is selected more than once (part S and part M)")]
  #[case::planned_output(  ("out/gtr.json",     "out/M.fasta"),  "Output destination 'out/gtr.json' is selected more than once (--output-gtr and part S)")]
  #[trace]
  fn test_output_plan_expansion_rejects_shared_destination(
    #[case] (path_s, path_m): (&str, &str),
    #[case] expected: &str,
  ) -> Result<(), Report> {
    let resolved = plan(&helpers::request(
      vec![OutputSelection::Gtr],
      btreemap! {
        OutputSelection::Gtr => PathBuf::from("out/gtr.json"),
        OutputSelection::ReconstructedAaFasta => PathBuf::from("out/aa.fasta"),
      },
    ))?;

    assert_error!(
      resolved.ensure_unique_with_expansion(
        OutputSelection::ReconstructedAaFasta,
        [(o!("part S"), Path::new(path_s)), (o!("part M"), Path::new(path_m))],
      ),
      expected
    );
    Ok(())
  }

  #[test]
  fn test_output_plan_expansion_replaces_the_expanded_output() -> Result<(), Report> {
    let resolved = plan(&helpers::request(
      vec![OutputSelection::Gtr],
      btreemap! { OutputSelection::ReconstructedAaFasta => PathBuf::from("out/aa.fasta") },
    ))?;

    resolved.ensure_unique_with_expansion(
      OutputSelection::ReconstructedAaFasta,
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

  mod helpers {
    use crate::output_plan::{CommandKind, OutputPlanRequest, OutputSelection};
    use std::collections::BTreeMap;
    use std::path::PathBuf;
    use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
    use treetime::progress::{LogLevel, LogSink};

    pub(super) fn request(
      selection: Vec<OutputSelection>,
      non_tree_overrides: BTreeMap<OutputSelection, PathBuf>,
    ) -> OutputPlanRequest {
      OutputPlanRequest {
        command: CommandKind::Ancestral,
        output_all: Some(PathBuf::from("out")),
        nwk_styles: vec![],
        selection,
        tree_overrides: BTreeMap::new(),
        non_tree_overrides,
      }
    }

    pub(super) fn homoplasy_request(non_tree_overrides: BTreeMap<OutputSelection, PathBuf>) -> OutputPlanRequest {
      OutputPlanRequest {
        command: CommandKind::Homoplasy,
        output_all: Some(PathBuf::from("out")),
        nwk_styles: vec![],
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
