#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::config::cli_flags::annotate_cli_flags;
  use crate::config::properties::{CLI_FLAG_KEY, PathRole, setting_properties};
  use clap::Arg;
  use helpers::{annotated_schema, path_hinted_args};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use schemars::Schema;
  use serde_json::{Value, json};
  use std::collections::BTreeSet;
  use strum::IntoEnumIterator;
  use treetime_utils::assert_error;

  #[test]
  fn test_schema_annotations_every_setting_of_every_command_has_a_flag() {
    for command in AppCommand::iter() {
      let schema = annotated_schema(command);
      for setting in setting_properties(&schema).unwrap() {
        let flag = schema.pointer(&setting.schema_pointer).unwrap()[CLI_FLAG_KEY].as_str();
        assert!(
          flag.is_some_and(|flag| flag.starts_with("--")),
          "{command}: `{}` has no flag",
          setting.key
        );
      }
    }
  }

  #[test]
  fn test_schema_annotations_every_setting_of_every_command_has_a_help_heading() {
    for command in AppCommand::iter() {
      let cli = command.cli_command();
      let missing = setting_properties(command.config_schema().as_value())
        .unwrap()
        .into_iter()
        .filter(|setting| {
          cli
            .get_arguments()
            .find(|arg| arg.get_id() == setting.key.as_str())
            .and_then(Arg::get_help_heading)
            .is_none()
        })
        .map(|setting| setting.key)
        .collect_vec();
      assert_eq!(
        Vec::<String>::new(),
        missing,
        "{command}: settings without a help heading"
      );
    }
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::timetree_input(    AppCommand::Timetree,  "tree",                        "Input data")]
  #[case::timetree_clock(    AppCommand::Timetree,  "clock_rate",                  "Molecular clock")]
  #[case::timetree_dating(   AppCommand::Timetree,  "max_iter",                    "Dating")]
  #[case::optimize_iter(     AppCommand::Optimize,  "max_iter",                    "Branch lengths")]
  #[case::prune_threshold(   AppCommand::Prune,     "prune_short",                 "Pruning")]
  #[case::clock_flattened(   AppCommand::Clock,     "reroot",                      "Rooting")]
  #[case::shared_alignment(  AppCommand::Ancestral, "alignment",                   "Input data")]
  #[case::shared_order(      AppCommand::Mugration, "topology_order",              "Tree ordering")]
  #[case::output_path(       AppCommand::Timetree,  "output_clock_model",          "Output")]
  #[case::hidden_plot(       AppCommand::Timetree,  "plot_tree",                   "Plots")]
  #[trace]
  fn test_schema_annotations_help_heading_of_setting(
    #[case] command: AppCommand,
    #[case] id: &str,
    #[case] expected: &str,
  ) {
    let cli = command.cli_command();
    let arg = cli.get_arguments().find(|arg| arg.get_id() == id).unwrap();
    assert_eq!(Some(expected), arg.get_help_heading());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::renamed_long(AppCommand::Timetree,  "metadata",  "--metadata")]
  #[case::flattened(   AppCommand::Ancestral, "alignment", "--alignment")]
  #[trace]
  fn test_schema_annotations_flag_of_setting(
    #[case] command: AppCommand,
    #[case] key: &str,
    #[case] expected: &str,
  ) {
    let schema = annotated_schema(command);
    let setting = setting_properties(&schema)
      .unwrap()
      .into_iter()
      .find(|setting| setting.key == key)
      .unwrap();
    assert_eq!(
      Some(expected),
      schema.pointer(&setting.schema_pointer).unwrap()[CLI_FLAG_KEY].as_str()
    );
  }

  #[test]
  fn test_schema_annotations_path_settings_match_the_path_arguments_of_the_cli() {
    for command in AppCommand::iter() {
      let schema = command.config_schema();
      let annotated: BTreeSet<String> = setting_properties(schema.as_value())
        .unwrap()
        .into_iter()
        .filter(|setting| setting.path_role.is_some())
        .map(|setting| setting.key)
        .collect();
      assert_eq!(path_hinted_args(&command.cli_command()), annotated, "{command}");
    }
  }

  #[test]
  fn test_schema_annotations_output_paths_are_the_output_and_plot_settings() {
    for command in AppCommand::iter() {
      for setting in setting_properties(command.config_schema().as_value()).unwrap() {
        let key = setting.key.as_str();
        let is_output_key = key.starts_with("output_") || key.starts_with("plot_");
        let is_output_role = setting.path_role == Some(PathRole::Output);
        assert_eq!(
          is_output_key && setting.path_role.is_some(),
          is_output_role,
          "{command}: `{key}` has role {:?}",
          setting.path_role
        );
      }
    }
  }

  #[test]
  fn test_schema_annotations_ancestral_path_roles() {
    let roles: BTreeSet<(String, String)> = setting_properties(AppCommand::Ancestral.config_schema().as_value())
      .unwrap()
      .into_iter()
      .filter_map(|setting| {
        let role = setting.path_role?;
        (role != PathRole::Output).then(|| (setting.key, format!("{role:?}")))
      })
      .collect();
    let expected: BTreeSet<(String, String)> = [
      ("aa_root_sequence", "Input"),
      ("alignment", "Input"),
      ("annotation", "Input"),
      ("custom_gtr", "Input"),
      ("topology_order_target_file", "Input"),
      ("translations", "InputTemplate"),
      ("tree", "Input"),
      ("vcf_reference", "Input"),
    ]
    .into_iter()
    .map(|(key, role)| (key.to_owned(), role.to_owned()))
    .collect();
    assert_eq!(expected, roles);
  }

  #[test]
  fn test_schema_annotations_setting_without_flag_is_an_error() {
    let mut schema = Schema::try_from(json!({
      "type": "object",
      "properties": { "not_a_flag": { "type": "string" } }
    }))
    .unwrap();
    let command = AppCommand::Prune.cli_command();
    assert_error!(
      annotate_cli_flags(&mut schema, &command),
      "config key `not_a_flag` has no command-line flag"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::definition(         json!({ "$ref": "#/$defs/Pair" }))]
  #[case::optional_definition(json!({ "anyOf": [{ "$ref": "#/$defs/Pair" }, { "type": "null" }] }))]
  #[case::inline_object(      json!({ "type": "object", "properties": { "x": { "type": "number" } } }))]
  #[trace]
  fn test_schema_annotations_setting_group_is_an_error(#[case] group: Value) {
    let schema = json!({
      "type": "object",
      "properties": { "a": group, "b": { "type": "number" } },
      "$defs": { "Pair": { "type": "object", "properties": { "x": { "type": "number" } } } }
    });
    assert_error!(
      setting_properties(&schema),
      "config key `a` holds a group of settings: every setting of a command must be a top-level key of its config"
    );
  }

  #[test]
  fn test_schema_annotations_settings_are_the_top_level_keys() {
    let schema = json!({
      "type": "object",
      "properties": {
        "$schema": { "type": "string" },
        "mode": { "$ref": "#/$defs/Mode" },
        "vars": { "type": "object", "additionalProperties": { "type": "string" } },
        "rate": { "anyOf": [{ "type": "number" }, { "type": "null" }] }
      },
      "$defs": { "Mode": { "type": "string", "enum": ["one", "two"] } }
    });
    assert_eq!(
      vec!["mode", "vars", "rate"],
      setting_properties(&schema)
        .unwrap()
        .into_iter()
        .map(|setting| setting.key)
        .collect_vec()
    );
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::config::cli_flags::annotate_cli_flags;
    use clap::{Command, ValueHint};
    use serde_json::Value;
    use std::collections::BTreeSet;

    pub(super) fn annotated_schema(command: AppCommand) -> Value {
      let mut schema = command.config_schema();
      annotate_cli_flags(&mut schema, &command.cli_command()).unwrap();
      schema.as_value().clone()
    }

    pub(super) fn path_hinted_args(command: &Command) -> BTreeSet<String> {
      command
        .get_arguments()
        .filter(|arg| {
          matches!(
            arg.get_value_hint(),
            ValueHint::FilePath | ValueHint::DirPath | ValueHint::AnyPath
          )
        })
        .map(|arg| arg.get_id().to_string())
        .filter(|id| id != "config")
        .collect()
    }
  }
}
