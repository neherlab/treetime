#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::command_config::CommandConfig;
  use crate::json_value::JsonValue;
  use crate::runs::setting_differences::{SettingDifference, setting_differences};
  use helpers::{input, record};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::path::PathBuf;
  use treetime_utils::o;

  #[test]
  fn test_setting_differences_of_identical_runs_are_empty() {
    let first = record(&json!({ "clock_rate": 0.001 }), vec![input("tree", "t.nwk", "aa")]);
    let second = record(&json!({ "clock_rate": 0.001 }), vec![input("tree", "t.nwk", "aa")]);
    assert_eq!(
      Vec::<SettingDifference>::new(),
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_lists_changed_values_in_catalog_order() {
    let first = record(&json!({ "clock_rate": 0.001, "max_iter": 2 }), vec![]);
    let second = record(&json!({ "clock_rate": 0.002, "max_iter": 5 }), vec![]);
    assert_eq!(
      vec![
        SettingDifference::Setting {
          key: o!("clock_rate"),
          first: Some(JsonValue(json!(0.001))),
          second: Some(JsonValue(json!(0.002))),
        },
        SettingDifference::Setting {
          key: o!("max_iter"),
          first: Some(JsonValue(json!(2))),
          second: Some(JsonValue(json!(5))),
        },
      ],
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_nested_setting_is_compared_by_its_key_path() {
    let mut first = record(&json!({}), vec![]);
    first.config = CommandConfig::from_settings(AppCommand::Clock, &json!({})).unwrap();
    let mut second = first.clone();
    second.config =
      CommandConfig::from_settings(AppCommand::Clock, &json!({ "branch_split": { "method": "brent" } })).unwrap();
    assert_eq!(
      vec![SettingDifference::Setting {
        key: o!("branch_split.method"),
        first: Some(JsonValue(json!("grid"))),
        second: Some(JsonValue(json!("brent"))),
      }],
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_moved_input_with_the_same_contents() {
    let first = record(&json!({}), vec![input("tree", "a/t.nwk", "aa")]);
    let second = record(&json!({}), vec![input("tree", "b/t.nwk", "aa")]);
    assert_eq!(
      vec![SettingDifference::Input {
        key: o!("tree"),
        first: vec![PathBuf::from("a/t.nwk")],
        second: vec![PathBuf::from("b/t.nwk")],
        same_content: true,
      }],
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_changed_input_contents_at_the_same_path() {
    let first = record(&json!({}), vec![input("metadata", "m.tsv", "aa")]);
    let second = record(&json!({}), vec![input("metadata", "m.tsv", "bb")]);
    assert_eq!(
      vec![SettingDifference::Input {
        key: o!("metadata"),
        first: vec![PathBuf::from("m.tsv")],
        second: vec![PathBuf::from("m.tsv")],
        same_content: false,
      }],
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_input_order_of_the_same_files_changes_paths_only() {
    let files = [input("alignment", "a.fasta", "aa"), input("alignment", "b.fasta", "bb")];
    let first = record(&json!({}), files.to_vec());
    let second = record(&json!({}), files.iter().rev().cloned().collect());
    assert_eq!(
      vec![SettingDifference::Input {
        key: o!("alignment"),
        first: vec![PathBuf::from("a.fasta"), PathBuf::from("b.fasta")],
        second: vec![PathBuf::from("b.fasta"), PathBuf::from("a.fasta")],
        same_content: true,
      }],
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_ignore_output_paths() {
    let first = record(&json!({ "output_all": "out", "output_tracelog": "a.csv" }), vec![]);
    let second = record(
      &json!({ "output_all": "elsewhere", "output_tracelog": "b.csv" }),
      vec![],
    );
    assert_eq!(
      Vec::<SettingDifference>::new(),
      setting_differences(&first, &second).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_serialize_with_kind_tag() {
    let difference = SettingDifference::Input {
      key: o!("tree"),
      first: vec![PathBuf::from("a.nwk")],
      second: vec![PathBuf::from("b.nwk")],
      same_content: false,
    };
    assert_eq!(
      json!({ "kind": "input", "key": "tree", "first": ["a.nwk"], "second": ["b.nwk"], "same_content": false }),
      serde_json::to_value(difference).unwrap()
    );
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::command_config::CommandConfig;
    use crate::job::JobId;
    use crate::runs::headline::RunHeadline;
    use crate::runs::record::{RunInput, RunRecord, RunStatus};
    use chrono::DateTime;
    use serde_json::Value;
    use std::path::PathBuf;

    pub(super) fn record(settings: &Value, inputs: Vec<RunInput>) -> RunRecord {
      RunRecord {
        id: JobId::random(),
        title: "run".to_owned(),
        config: CommandConfig::from_settings(AppCommand::Timetree, settings).unwrap(),
        status: RunStatus::Ok,
        pinned: false,
        created_at: DateTime::UNIX_EPOCH,
        started_at: None,
        finished_at: None,
        duration_seconds: None,
        treetime_version: "1.0.0".to_owned(),
        inputs,
        config_hash: None,
        changed_settings: vec![],
        headline: RunHeadline::default(),
        output_files: vec![],
        error: None,
      }
    }

    pub(super) fn input(setting: &str, path: &str, sha256: &str) -> RunInput {
      RunInput {
        setting: setting.to_owned(),
        path: PathBuf::from(path),
        size: 1,
        sha256: sha256.to_owned(),
      }
    }
  }
}
