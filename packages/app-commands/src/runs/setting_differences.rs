use crate::config::catalog::{SettingRole, SettingSpec, command_settings};
use crate::config::settings::setting_ref;
use crate::runs::record::{RunInput, RunRecord};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::PathBuf;

pub fn setting_differences(first: &RunRecord, second: &RunRecord) -> Result<Vec<SettingDifference>, Report> {
  Ok(
    command_settings(first.command)?
      .settings
      .iter()
      .filter_map(|spec| match spec.role {
        SettingRole::Output => None,
        SettingRole::Input | SettingRole::InputTemplate => input_difference(spec, first, second),
        SettingRole::Setting => {
          let first_value = setting_value(first, spec);
          let second_value = setting_value(second, spec);
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "kebab-case")]
pub enum SettingDifference {
  /// A value that controls the analysis.
  Setting {
    /// Key path of the setting joined with `.`.
    key: String,
    /// Value in the first run.
    first: Value,
    /// Value in the second run.
    second: Value,
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

fn setting_value(record: &RunRecord, spec: &SettingSpec) -> Value {
  setting_ref(&record.config, &spec.path)
    .cloned()
    .unwrap_or_else(|| spec.default_value.clone())
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
