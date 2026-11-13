#[cfg(test)]
mod tests {
  use crate::config::catalog::{SettingPurpose, command_settings};
  use crate::runs::setting_differences::{ComparedConfig, DifferenceValues, SettingDifference, setting_differences};
  use app_commands::command::AppCommand;
  use app_commands::json_value::JsonValue;
  use helpers::{differences, input, unread};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::path::PathBuf;
  use treetime_utils::io::json::to_json_value;
  use treetime_utils::o;

  #[test]
  fn test_setting_differences_of_identical_configurations_are_empty() {
    assert_eq!(
      Vec::<SettingDifference>::new(),
      differences(
        (&json!({ "clock_rate": 0.001 }), &[input("tree", "t.nwk", "aa")]),
        (&json!({ "clock_rate": 0.001 }), &[input("tree", "t.nwk", "aa")]),
      )
    );
  }

  #[test]
  fn test_setting_differences_list_changed_values_with_their_labels_in_catalog_order() {
    assert_eq!(
      vec![
        SettingDifference {
          key: o!("clock_rate"),
          label: o!("Clock rate"),
          group: o!("Molecular clock"),
          purpose: SettingPurpose::Analysis,
          values: DifferenceValues::Setting {
            first: None,
            second: Some(JsonValue(json!(0.002))),
          },
        },
        SettingDifference {
          key: o!("max_iter"),
          label: o!("Iterations"),
          group: o!("Branch lengths"),
          purpose: SettingPurpose::Analysis,
          values: DifferenceValues::Setting {
            first: Some(JsonValue(json!(2))),
            second: Some(JsonValue(json!(5))),
          },
        },
      ],
      differences(
        (&json!({ "max_iter": 2 }), &[]),
        (&json!({ "clock_rate": 0.002, "max_iter": 5 }), &[]),
      )
    );
  }

  #[test]
  fn test_setting_differences_compare_inputs_by_content_and_not_by_path() {
    let files = [input("alignment", "a.fasta", "aa"), input("alignment", "b.fasta", "bb")];
    let reordered = [
      input("alignment", "elsewhere/b.fasta", "bb"),
      input("alignment", "inputs/a.fasta", "aa"),
    ];
    assert_eq!(
      Vec::<SettingDifference>::new(),
      differences((&json!({}), &files), (&json!({}), &reordered))
    );
  }

  #[test]
  fn test_setting_differences_list_inputs_with_changed_contents() {
    assert_eq!(
      vec![SettingDifference {
        key: o!("metadata"),
        label: o!("Metadata"),
        group: o!("Input data"),
        purpose: SettingPurpose::Data,
        values: DifferenceValues::Input {
          first: vec![PathBuf::from("m.tsv")],
          second: vec![PathBuf::from("m.tsv")],
        },
      }],
      differences(
        (&json!({}), &[input("metadata", "m.tsv", "aa")]),
        (&json!({}), &[input("metadata", "m.tsv", "bb")]),
      )
    );
  }

  #[test]
  fn test_setting_differences_count_an_unreadable_input_as_changed() {
    assert_eq!(
      vec![o!("tree")],
      differences(
        (&json!({}), &[input("tree", "t.nwk", "aa")]),
        (&json!({}), &[unread("tree", "t.nwk")]),
      )
      .into_iter()
      .map(|difference| difference.key)
      .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_setting_differences_ignore_the_settings_that_the_app_sets_for_every_run() {
    assert_eq!(
      Vec::<SettingDifference>::new(),
      differences(
        (
          &json!({ "output_all": "out", "output_tracelog": "a.csv", "output_selection": ["auspice", "tracelog"] }),
          &[],
        ),
        (&json!({ "output_all": "elsewhere", "output_tracelog": "b.csv" }), &[]),
      )
    );
  }

  #[test]
  fn test_setting_differences_serialize_flat_with_a_kind_tag() {
    let difference = SettingDifference {
      key: o!("tree"),
      label: o!("Tree"),
      group: o!("Input data"),
      purpose: SettingPurpose::Data,
      values: DifferenceValues::Input {
        first: vec![PathBuf::from("a.nwk")],
        second: vec![PathBuf::from("b.nwk")],
      },
    };
    assert_eq!(
      json!({
        "key": "tree",
        "label": "Tree",
        "group": "Input data",
        "purpose": "data",
        "kind": "input",
        "first": ["a.nwk"],
        "second": ["b.nwk"],
      }),
      to_json_value(&difference).unwrap()
    );
  }

  #[test]
  fn test_setting_differences_of_an_unset_setting_omit_its_value() {
    let catalog = command_settings(AppCommand::Timetree).unwrap();
    let first = json!({ "coalescent": "skyline" });
    let second = json!({ "coalescent": null });
    assert_eq!(
      json!([{
        "key": "coalescent",
        "label": "Coalescent prior",
        "group": "Coalescent prior",
        "purpose": "analysis",
        "kind": "setting",
        "first": "skyline",
      }]),
      to_json_value(&setting_differences(
        &catalog,
        &ComparedConfig {
          settings: first.as_object().unwrap(),
          inputs: &[],
        },
        &ComparedConfig {
          settings: second.as_object().unwrap(),
          inputs: &[],
        },
      ))
      .unwrap()
    );
  }

  mod helpers {
    use crate::command::config_from_settings;
    use crate::config::catalog::command_settings;
    use crate::runs::setting_differences::{ComparedConfig, ComparedInput, SettingDifference, setting_differences};
    use app_commands::command::AppCommand;
    use serde_json::Value;
    use std::path::PathBuf;

    pub(super) fn differences(
      (first, first_inputs): (&Value, &[ComparedInput]),
      (second, second_inputs): (&Value, &[ComparedInput]),
    ) -> Vec<SettingDifference> {
      let catalog = command_settings(AppCommand::Timetree).unwrap();
      let resolved = |settings: &Value| {
        config_from_settings(AppCommand::Timetree, settings)
          .unwrap()
          .settings()
          .unwrap()
      };
      setting_differences(
        &catalog,
        &ComparedConfig {
          settings: &resolved(first),
          inputs: first_inputs,
        },
        &ComparedConfig {
          settings: &resolved(second),
          inputs: second_inputs,
        },
      )
    }

    pub(super) fn input(setting: &str, path: &str, sha256: &str) -> ComparedInput {
      ComparedInput {
        setting: setting.to_owned(),
        path: PathBuf::from(path),
        sha256: Some(sha256.to_owned()),
      }
    }

    pub(super) fn unread(setting: &str, path: &str) -> ComparedInput {
      ComparedInput {
        setting: setting.to_owned(),
        path: PathBuf::from(path),
        sha256: None,
      }
    }
  }
}
