use crate::command::AppCommandExt;
use crate::config::catalog::{SettingRole, SettingSpec};
use app_commands::command::AppCommand;
use app_commands::commands::shared::reroot::RerootMethodCli;
use deser::Serialize;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use treetime_utils::io::json::to_json_value;

const AUTOMATIC: [(&str, Automatic); 5] = [
  ("clock_rate", Automatic::Text("Estimated from the data")),
  ("clock_filter", Automatic::ClockFilter),
  (
    "coalescent",
    Automatic::Text("None: the node times have no coalescent prior"),
  ),
  ("reroot", Automatic::Reroot),
  (
    "seed",
    Automatic::Text("Random: drawn for the run and written to the log"),
  ),
];

pub fn effective_values(
  command: AppCommand,
  settings: &[SettingSpec],
  config: &Map<String, Value>,
) -> Result<Vec<EffectiveValue>, Report> {
  settings
    .iter()
    .filter(|spec| spec.role == SettingRole::Setting && is_unset(spec, config))
    .map(|spec| {
      let text = match AUTOMATIC.iter().find(|(key, _)| *key == spec.key) {
        Some((_, automatic)) => automatic.text(command, spec, config)?,
        None => value_text(
          spec
            .default_value
            .as_ref()
            .or(spec.documented_default.as_ref())
            .map(|value| &value.0),
        ),
      };
      Ok(EffectiveValue {
        key: spec.key.clone(),
        text,
      })
    })
    .try_collect()
}

/// The value that a run uses for a setting that the configuration leaves at its default.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
pub struct EffectiveValue {
  /// Key of the setting.
  pub key: String,
  /// The value as a short text, for example `Automatic: least-squares rerooting`.
  pub text: String,
}

fn is_unset(spec: &SettingSpec, config: &Map<String, Value>) -> bool {
  let value = config.get(&spec.key).filter(|value| !value.is_null());
  let default = spec
    .default_value
    .as_ref()
    .map(|value| &value.0)
    .filter(|value| !value.is_null());
  value == default
}

enum Automatic {
  Text(&'static str),
  ClockFilter,
  Reroot,
}

impl Automatic {
  fn text(&self, command: AppCommand, spec: &SettingSpec, config: &Map<String, Value>) -> Result<String, Report> {
    match self {
      Self::Text(text) => Ok((*text).to_owned()),
      Self::ClockFilter => Ok(clock_filter_text(
        spec.default_value.as_ref().and_then(|value| value.0.as_f64()),
      )),
      Self::Reroot => reroot_text(command, config),
    }
  }
}

fn clock_filter_text(threshold: Option<f64>) -> String {
  match threshold {
    Some(threshold) if threshold > 0.0 => format!("{threshold} interquartile ranges"),
    _ => "Off: all samples stay in the regression".to_owned(),
  }
}

fn reroot_text(command: AppCommand, config: &Map<String, Value>) -> Result<String, Report> {
  let keep_root = config.get("keep_root").and_then(Value::as_bool).unwrap_or(false);
  let outgroup = config
    .get("outgroup")
    .and_then(Value::as_array)
    .is_some_and(|tips| !tips.is_empty());
  if keep_root {
    return Ok("Not used: the input root is kept".to_owned());
  }
  if outgroup {
    return Ok("Not used: the tree is rooted on the outgroup".to_owned());
  }
  if !command.reroots_by_default() {
    return Ok("Not used: the input root is kept".to_owned());
  }
  let method = to_json_value(&RerootMethodCli::default())?;
  Ok(format!("Automatic: {} rerooting", value_text(Some(&method))))
}

pub(crate) fn value_text(value: Option<&Value>) -> String {
  match value {
    None | Some(Value::Null) => "Not set".to_owned(),
    Some(Value::Bool(true)) => "On".to_owned(),
    Some(Value::Bool(false)) => "Off".to_owned(),
    Some(Value::Number(number)) => number
      .as_i64()
      .map(|integer| integer.to_string())
      .or_else(|| number.as_f64().map(|float| float.to_string()))
      .unwrap_or_else(|| number.to_string()),
    Some(Value::String(text)) => text.clone(),
    Some(Value::Array(items)) if items.is_empty() => "None".to_owned(),
    Some(Value::Array(items)) => items.iter().map(|item| value_text(Some(item))).join(", "),
    Some(object @ Value::Object(_)) => object.to_string(),
  }
}
