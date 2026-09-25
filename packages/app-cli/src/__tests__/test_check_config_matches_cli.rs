#[cfg(test)]
mod tests {
  use app_commands::command::AppCommand;
  use helpers::{check_config_error, cli_error};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown_field(    AppCommand::Ancestral, "tree: t.nwk\ndefinitely_not_a_real_field: 1\n")]
  #[case::unknown_nested(   AppCommand::Clock,     "tree: t.nwk\nmetadata: m.tsv\nbranch_split:\n  n_pointz: 3\n")]
  #[case::wrong_enum(       AppCommand::Ancestral, "tree: t.nwk\nmethod_anc: margnal\n")]
  #[case::wrong_type(       AppCommand::Timetree,  "clock_rate: fast\n")]
  #[case::duplicate_key(    AppCommand::Prune,     "tree: a.nwk\ntree: b.nwk\n")]
  #[case::missing_required( AppCommand::Mugration, "tree: t.nwk\n")]
  #[case::missing_tree(     AppCommand::Optimize,  "seed: 3\n")]
  #[trace]
  fn test_check_config_matches_cli_error(#[case] command: AppCommand, #[case] text: &str) {
    assert_eq!(cli_error(command, text), check_config_error(command, text));
  }

  mod helpers {
    use crate::cli::treetime_cli::treetime_parse_cli_args;
    use crate::run::run_command;
    use app_commands::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
    use app_commands::command::AppCommand;
    use serde_json::Map;
    use std::fs;
    use tempfile::tempdir;
    use treetime::progress::NoopProgress;
    use treetime_utils::error::report_to_string;

    pub(super) fn cli_error(command: AppCommand, text: &str) -> String {
      let dir = tempdir().unwrap();
      let path = dir.path().join("config.yaml");
      fs::write(&path, text).unwrap();
      let argv = [
        "treetime".to_owned(),
        command.to_string(),
        "--config".to_owned(),
        path.to_string_lossy().into_owned(),
      ];
      let result = treetime_parse_cli_args(argv).and_then(|args| run_command(args.command, &NoopProgress));
      match result {
        Ok(()) => panic!("the CLI accepted the config"),
        Err(report) => report_to_string(&report),
      }
    }

    pub(super) fn check_config_error(command: AppCommand, text: &str) -> String {
      let response = check_config(&CheckConfigRequest {
        command,
        text: text.to_owned(),
        inputs: Map::new(),
        input_facts: None,
      });
      match response {
        CheckConfigResponse::Valid { .. } => panic!("check-config accepted the config"),
        CheckConfigResponse::Invalid { message, causes, .. } => [vec![message], causes].concat().join(": "),
      }
    }
  }
}
