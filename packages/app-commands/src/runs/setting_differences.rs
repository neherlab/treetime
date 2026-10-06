use crate::config::catalog::{SettingRole, SettingSpec, command_settings};
use crate::config::settings::setting_ref;
use crate::json_value::JsonValue;
use crate::runs::record::{RunInput, RunRecord};
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::path::PathBuf;
use treetime_schema::skip_serializing_optionals;

pub fn setting_differences(first: &RunRecord, second: &RunRecord) -> Result<Vec<SettingDifference>, Report> {
  let first_settings = first.config.settings()?;
  let second_settings = second.config.settings()?;
  Ok(
    command_settings(first.config.command())?
      .settings
      .iter()
      .filter_map(|spec| match spec.role {
        SettingRole::Output => None,
        SettingRole::Input | SettingRole::InputTemplate => input_difference(spec, first, second),
        SettingRole::Setting => {
          let first_value = setting_value(&first_settings, spec);
          let second_value = setting_value(&second_settings, spec);
          (first_value != second_value).then(|| SettingDifference::Setting {
            key: spec.key.clone(),
            first: first_value,
            second: second_value,
          })
        },
      })
      .collect(),
  )
}

/// A setting whose value differs between two runs.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum SettingDifference {
  /// A value that controls the analysis.
  Setting {
    /// Key path of the setting joined with `.`.
    key: String,
    /// Value in the first run; absent when the first run does not set the setting.
    first: Option<JsonValue>,
    /// Value in the second run; absent when the second run does not set the setting.
    second: Option<JsonValue>,
  },
  /// Input files, which differ in their paths, their contents, or both.
  Input {
    /// Key path of the setting that names the files.
    key: String,
    /// Paths of the files the first run read.
    first: Vec<PathBuf>,
    /// Paths of the files the second run read.
    second: Vec<PathBuf>,
    /// Whether the files of both runs have the same contents.
    same_content: bool,
  },
}

fn setting_value(settings: &Map<String, Value>, spec: &SettingSpec) -> Option<JsonValue> {
  setting_ref(settings, &spec.path).cloned().map(JsonValue)
}

fn input_difference(spec: &SettingSpec, first: &RunRecord, second: &RunRecord) -> Option<SettingDifference> {
  let first_inputs = inputs_of(first, &spec.key);
  let second_inputs = inputs_of(second, &spec.key);
  let first_paths = first_inputs.iter().map(|input| input.path.clone()).collect_vec();
  let second_paths = second_inputs.iter().map(|input| input.path.clone()).collect_vec();
  let same_content = contents(&first_inputs) == contents(&second_inputs);
  (first_paths != second_paths || !same_content).then(|| SettingDifference::Input {
    key: spec.key.clone(),
    first: first_paths,
    second: second_paths,
    same_content,
  })
}

fn inputs_of<'a>(record: &'a RunRecord, key: &str) -> Vec<&'a RunInput> {
  record.inputs.iter().filter(|input| input.setting == key).collect()
}

fn contents(inputs: &[&RunInput]) -> Vec<String> {
  inputs.iter().map(|input| input.sha256.clone()).sorted().collect()
}
