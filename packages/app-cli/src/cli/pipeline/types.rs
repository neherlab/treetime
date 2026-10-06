use app_commands::command::AppCommand;
use app_commands::command_config::CommandConfig;
use app_commands::config::schema::SCHEMA_KEY;
use app_commands::config::suggest::suggestion_suffix;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde_json::{Map, Value, json};
use std::borrow::Cow;
use std::str::FromStr;
use strum::{IntoEnumIterator, VariantNames};
use treetime_utils::{make_error, make_report};

/// The whole pipeline in typed form, used for schema generation.
///
/// The loader does not deserialize into this type directly: `vars`, `output_all`, and each step's
/// outputs are resolved in a staged, backward-only pass (interpolation depends on earlier results),
/// after which the typed steps are assembled. This type fixes the on-disk shape that the staged pass
/// and the generated schema must agree on.
#[derive(Debug, JsonSchema)]
#[serde(deny_unknown_fields)]
#[expect(
  dead_code,
  reason = "the type only describes the file shape for the generated schema; the loader reads files in a staged pass"
)]
pub(crate) struct Pipeline {
  #[serde(rename = "$schema", skip_serializing_if = "Option::is_none")]
  schema_ref: Option<String>,

  #[serde(default, skip_serializing_if = "Map::is_empty")]
  vars: Map<String, Value>,

  #[serde(skip_serializing_if = "Option::is_none")]
  output_all: Option<String>,

  steps: Vec<PipelineStep>,
}

/// One named step: a stable id plus exactly one command invocation.
///
/// The `name` is explicit (not the command name) because `--steps=` selection and
/// `{{ steps.<name>... }}` references need stable ids and must allow the same command twice.
#[derive(Debug, JsonSchema)]
#[expect(
  dead_code,
  reason = "the type only describes the file shape for the generated schema; the loader reads steps in a staged pass"
)]
pub(crate) struct PipelineStep {
  name: String,
  #[serde(flatten)]
  command: PipelineStepCommand,
}

/// A single analysis command invocation within a pipeline: the command name as the key, and the command's settings,
/// in the shape the per-command `--config` accepts, as the value.
#[derive(Debug)]
pub(crate) struct PipelineStepCommand;

impl JsonSchema for PipelineStepCommand {
  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("PipelineStepCommand")
  }

  fn json_schema(generator: &mut SchemaGenerator) -> Schema {
    let variants = AppCommand::iter()
      .map(|command| {
        let name: &str = command.into();
        json!({
          "type": "object",
          "properties": { name: command.config_subschema(generator) },
          "required": [name],
        })
      })
      .collect_vec();
    json_schema!({ "oneOf": variants })
  }
}

pub(crate) fn step_config(tag: &str, payload: &Value) -> Result<CommandConfig, Report> {
  let command = AppCommand::from_str(tag).map_err(|_unknown| {
    make_report!(
      "unknown command `{tag}`; {}",
      suggestion_suffix(tag, AppCommand::VARIANTS)
    )
  })?;
  command
    .command_config(payload)
    .wrap_err_with(|| format!("When reading the settings of a `{command}` step"))
}

pub(crate) struct RawStep {
  pub name: String,
  pub tag: String,
  pub payload: Value,
}

impl RawStep {
  pub(crate) fn from_value(value: Value) -> Result<Self, Report> {
    let Value::Object(mut map) = value else {
      return make_error!("a pipeline step must be a mapping with a `name` and one command");
    };
    map.remove(SCHEMA_KEY);

    let name = match map.remove("name") {
      Some(Value::String(name)) => name,
      Some(_) => return make_error!("pipeline step `name` must be a string"),
      None => return make_error!("pipeline step is missing a `name`"),
    };

    let tags: Vec<String> = map.keys().cloned().collect();
    let (tag, payload) = match tags.as_slice() {
      [] => {
        return make_error!(
          "pipeline step `{name}` has no command; expected one of {}",
          commands_list()
        );
      },
      [tag] => (tag.clone(), map.remove(tag).unwrap_or(Value::Null)),
      _ => {
        return make_error!(
          "pipeline step `{name}` has more than one command ({}); a step runs exactly one command",
          tags.iter().sorted().map(|tag| format!("`{tag}`")).join(", ")
        );
      },
    };

    Ok(Self { name, tag, payload })
  }
}

pub(crate) fn commands_list() -> String {
  AppCommand::VARIANTS.iter().map(|tag| format!("`{tag}`")).join(", ")
}

#[cfg(test)]
mod tests {

  use helpers::{method_anc, parse_step};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_types_step_parses_ancestral_with_snake_case_fields() {
    let value = json!({
      "name": "anc",
      "ancestral": { "tree": "t.nwk", "method_anc": "marginal", "dense": true }
    });
    let step = parse_step(value).unwrap();
    assert_eq!(
      (o!("ancestral"), Some(o!("Marginal"))),
      (step.command().to_string(), method_anc(&step))
    );
  }

  #[test]
  fn test_types_step_parses_timetree_tag() {
    let step = parse_step(json!({ "name": "tt", "timetree": { "clock_rate": 0.003 } })).unwrap();
    assert_eq!("timetree", step.command().to_string());
  }

  #[test]
  fn test_types_step_rejects_utility_command_tag() {
    let result = parse_step(json!({ "name": "x", "debug": {} }));
    assert_error!(
      result,
      "in pipeline step `x`: unknown command `debug`; valid values: `ancestral`, `clock`, `homoplasy`, `mugration`, `optimize`, `prune`, `timetree`"
    );
  }

  #[test]
  fn test_types_step_suggests_closest_command_for_typo() {
    let result = parse_step(json!({ "name": "x", "timtree": {} }));
    assert_error!(
      result,
      "in pipeline step `x`: unknown command `timtree`; did you mean `timetree`? Valid values: `ancestral`, `clock`, `homoplasy`, `mugration`, `optimize`, `prune`, `timetree`"
    );
  }

  #[test]
  fn test_types_step_rejects_multiple_commands() {
    let result = parse_step(json!({ "name": "x", "timetree": {}, "clock": {} }));
    assert_error!(
      result,
      "pipeline step `x` has more than one command (`clock`, `timetree`); a step runs exactly one command"
    );
  }

  #[test]
  fn test_types_step_rejects_missing_command() {
    let result = parse_step(json!({ "name": "x" }));
    assert_error!(
      result,
      "pipeline step `x` has no command; expected one of `timetree`, `clock`, `ancestral`, `homoplasy`, `mugration`, `optimize`, `prune`"
    );
  }

  #[test]
  fn test_types_step_rejects_missing_name() {
    let result = parse_step(json!({ "timetree": {} }));
    assert_error!(result, "pipeline step is missing a `name`");
  }

  #[test]
  fn test_types_step_ignores_schema_key() {
    let value = json!({ "$schema": "./input-config-pipeline.schema.json", "name": "tt", "timetree": {} });
    let step = parse_step(value).unwrap();
    assert_eq!("timetree", step.command().to_string());
  }

  mod helpers {
    use super::super::*;
    use treetime_utils::make_report;

    pub(super) fn parse_step(value: Value) -> Result<CommandConfig, Report> {
      let RawStep { name, tag, payload } = RawStep::from_value(value)?;
      step_config(&tag, &payload).map_err(|err| make_report!("in pipeline step `{name}`: {err}"))
    }

    pub(super) fn method_anc(command: &CommandConfig) -> Option<String> {
      match command {
        CommandConfig::Ancestral(args) => Some(format!("{:?}", args.method_anc)),
        _ => None,
      }
    }
  }
}
