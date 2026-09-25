use crate::command::AppCommand;
use crate::config::properties::{
  CLI_FLAG_KEY, CLI_NUM_ARGS_KEY, CLI_VALUE_DELIMITER_KEY, CLI_VALUES_KEY, leaf_properties,
};
use crate::config::source::escape_pointer;
use clap::{Arg, Command};
use eyre::Report;
use itertools::Itertools;
use schemars::Schema;
use serde_json::{Map, Value, json};
use treetime_utils::{make_error, make_report};

pub fn annotated_config_schema(command: AppCommand) -> Result<Schema, Report> {
  let mut schema = command.config_schema();
  annotate_cli_flags(&mut schema, &command.cli_command())?;
  Ok(schema)
}

pub fn annotate_cli_flags(schema: &mut Schema, command: &Command) -> Result<(), Report> {
  let mut command = command.clone();
  command.build();
  let value = schema.as_value().clone();
  let leaves = leaf_properties(&value)?;
  let object = schema.ensure_object();
  let mut annotated = Value::Object(object.clone());
  for leaf in &leaves {
    let key_path = leaf.key_path.join(".");
    let arg = find_arg(&command, leaf.key())
      .filter(|arg| arg.get_long().is_some())
      .ok_or_else(|| make_report!("config key `{key_path}` has no command-line flag"))?;
    let property = value
      .pointer(&leaf.schema_pointer)
      .ok_or_else(|| make_report!("schema has no property at `{}`", leaf.schema_pointer))?;
    let annotations = cli_annotations(&value, property, arg, &key_path)?;
    let Some(Value::Object(property)) = annotated.pointer_mut(&leaf.schema_pointer) else {
      return make_error!("schema has no property at `{}`", leaf.schema_pointer);
    };
    property.extend(annotations);
  }
  let Value::Object(annotated) = annotated else {
    return make_error!("a command schema must be a JSON object");
  };
  *object = annotated;
  Ok(())
}

fn find_arg<'a>(command: &'a Command, id: &str) -> Option<&'a Arg> {
  command.get_arguments().find(|arg| arg.get_id() == id)
}

fn cli_annotations(schema: &Value, property: &Value, arg: &Arg, key_path: &str) -> Result<Map<String, Value>, Report> {
  let mut annotations = Map::new();
  if let Some(long) = arg.get_long() {
    annotations.insert(CLI_FLAG_KEY.to_owned(), Value::String(format!("--{long}")));
  }

  let (min, max) = if arg.get_action().takes_values() {
    arg.get_num_args().map_or((1, Some(1)), |range| {
      let max = range.max_values();
      (range.min_values(), (max != usize::MAX).then_some(max))
    })
  } else {
    (0, Some(0))
  };
  annotations.insert(CLI_NUM_ARGS_KEY.to_owned(), json!([min, max]));

  if let Some(delimiter) = arg.get_value_delimiter() {
    annotations.insert(CLI_VALUE_DELIMITER_KEY.to_owned(), Value::String(delimiter.to_string()));
  }

  let config_values = enum_values(schema, property);
  let cli_values = arg
    .get_possible_values()
    .iter()
    .map(|value| value.get_name().to_owned())
    .collect_vec();
  if !config_values.is_empty() && !cli_values.is_empty() {
    let mapping: Map<String, Value> = config_values
      .iter()
      .map(|config_value| {
        let cli_value = cli_values
          .iter()
          .find(|cli_value| spelling_key(cli_value) == spelling_key(config_value))
          .ok_or_else(|| {
            make_report!(
              "value `{config_value}` of config key `{key_path}` has no command-line spelling among {}",
              cli_values.join(", ")
            )
          })?;
        Ok((config_value.clone(), Value::String(cli_value.clone())))
      })
      .collect::<Result<_, Report>>()?;
    annotations.insert(CLI_VALUES_KEY.to_owned(), Value::Object(mapping));
  }

  Ok(annotations)
}

fn enum_values(schema: &Value, property: &Value) -> Vec<String> {
  let mut values = Vec::new();
  collect_enum_values(schema, property, &mut values);
  values.into_iter().unique().collect()
}

fn collect_enum_values(schema: &Value, node: &Value, values: &mut Vec<String>) {
  if let Some(reference) = node.get("$ref").and_then(Value::as_str)
    && let Some(def) = reference
      .strip_prefix("#/$defs/")
      .and_then(|name| schema.pointer(&format!("/$defs/{}", escape_pointer(name))))
  {
    collect_enum_values(schema, def, values);
  }
  if let Some(items) = node.get("enum").and_then(Value::as_array) {
    values.extend(items.iter().filter_map(Value::as_str).map(str::to_owned));
  }
  if let Some(constant) = node.get("const").and_then(Value::as_str) {
    values.push(constant.to_owned());
  }
  for key in ["oneOf", "anyOf", "allOf"] {
    for alternative in node.get(key).and_then(Value::as_array).into_iter().flatten() {
      collect_enum_values(schema, alternative, values);
    }
  }
  if let Some(items) = node.get("items") {
    collect_enum_values(schema, items, values);
  }
}

fn spelling_key(value: &str) -> String {
  value
    .chars()
    .filter(|c| *c != '-' && *c != '_')
    .flat_map(char::to_lowercase)
    .collect()
}
