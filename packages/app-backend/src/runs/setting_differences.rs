use crate::config::catalog::{CommandSettings, SettingPurpose, SettingRole, SettingSpec};
use crate::runs::record::RunInput;
use app_commands::json_value::JsonValue;
use deser::Serialize;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::path::PathBuf;
use treetime_schema::skip_serializing_optionals;

pub const OUTPUT_SELECTION_KEY: &str = "output_selection";

pub fn setting_differences(
  catalog: &CommandSettings,
  first: &ComparedConfig,
  second: &ComparedConfig,
) -> Vec<SettingDifference> {
  compared_settings(catalog)
    .filter_map(|spec| {
      let values = match spec.role {
        SettingRole::Output => None,
        SettingRole::Input | SettingRole::InputTemplate => input_values(spec, first.inputs, second.inputs),
        SettingRole::Setting => setting_values(spec, first.settings, second.settings),
      }?;
      Some(SettingDifference {
        key: spec.key.clone(),
        label: spec.label.clone(),
        group: spec.group.clone(),
        purpose: spec.purpose,
        values,
      })
    })
    .collect()
}

pub fn compared_settings(catalog: &CommandSettings) -> impl Iterator<Item = &SettingSpec> {
  catalog.settings.iter().filter(|spec| !is_set_by_app(spec))
}

pub fn compared_inputs(inputs: &[RunInput]) -> Vec<ComparedInput> {
  inputs
    .iter()
    .map(|input| ComparedInput {
      setting: input.setting.clone(),
      path: input.path.clone(),
      sha256: Some(input.sha256.clone()),
    })
    .collect()
}

/// A setting or an input whose resolved value differs between two configurations: a draft and a run, or two runs.
/// Inputs differ when their file contents differ; their paths do not count.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct SettingDifference {
  /// Key path of the setting joined with `.`.
  pub key: String,
  /// Name of the setting, as the setting catalog labels it, for example `Coalescent prior`.
  pub label: String,
  /// Group of the setting, as the headings of `treetime <command> --help` name it, for example `Coalescent prior`.
  pub group: String,
  /// What the setting is for.
  pub purpose: SettingPurpose,
  /// The values in the first and in the second configuration.
  #[schemars(flatten)]
  #[deser(flatten)]
  pub values: DifferenceValues,
}

/// Values of a setting that differs, in the first and in the second configuration.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(tag = "kind", rename_all = "kebab-case")]
#[deser(tag = "kind", rename_all = "kebab-case")]
pub enum DifferenceValues {
  /// A value that controls the analysis.
  Setting {
    /// Value in the first configuration; absent when it leaves the setting unset.
    first: Option<JsonValue>,
    /// Value in the second configuration; absent when it leaves the setting unset.
    second: Option<JsonValue>,
  },
  /// Input files with different contents.
  Input {
    /// Paths of the files of the first configuration.
    first: Vec<PathBuf>,
    /// Paths of the files of the second configuration.
    second: Vec<PathBuf>,
  },
}

pub struct ComparedConfig<'a> {
  pub settings: &'a Map<String, Value>,
  pub inputs: &'a [ComparedInput],
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ComparedInput {
  pub setting: String,
  pub path: PathBuf,
  pub sha256: Option<String>,
}

fn is_set_by_app(spec: &SettingSpec) -> bool {
  spec.role == SettingRole::Output || spec.key == OUTPUT_SELECTION_KEY
}

fn setting_values(
  spec: &SettingSpec,
  first: &Map<String, Value>,
  second: &Map<String, Value>,
) -> Option<DifferenceValues> {
  let first = setting_value(first, spec);
  let second = setting_value(second, spec);
  (first != second).then_some(DifferenceValues::Setting { first, second })
}

fn setting_value(settings: &Map<String, Value>, spec: &SettingSpec) -> Option<JsonValue> {
  settings
    .get(&spec.key)
    .filter(|value| !value.is_null())
    .cloned()
    .map(JsonValue)
}

fn input_values(spec: &SettingSpec, first: &[ComparedInput], second: &[ComparedInput]) -> Option<DifferenceValues> {
  let first = inputs_of(first, &spec.key);
  let second = inputs_of(second, &spec.key);
  let same_contents = matches!(
    (contents(&first), contents(&second)),
    (Some(first), Some(second)) if first == second
  );
  (!same_contents).then(|| DifferenceValues::Input {
    first: paths(&first),
    second: paths(&second),
  })
}

fn inputs_of<'a>(inputs: &'a [ComparedInput], key: &str) -> Vec<&'a ComparedInput> {
  inputs.iter().filter(|input| input.setting == key).collect()
}

fn contents(inputs: &[&ComparedInput]) -> Option<Vec<String>> {
  inputs
    .iter()
    .map(|input| input.sha256.clone())
    .collect::<Option<Vec<_>>>()
    .map(|contents| contents.into_iter().sorted().collect())
}

fn paths(inputs: &[&ComparedInput]) -> Vec<PathBuf> {
  inputs.iter().map(|input| input.path.clone()).collect()
}
