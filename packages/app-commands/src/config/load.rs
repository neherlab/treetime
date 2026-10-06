use crate::config::resolve_paths::resolve_config_paths;
use crate::config::schema::{SCHEMA_KEY, command_schema};
use crate::config::schema_check::schema_diagnostics;
use crate::config::source::{ConfigSource, parse_config_document, render_and_bail};
use deser::Serialize;
use eyre::Report;
use schemars::{JsonSchema, Schema};
use serde_json::Value;
use std::path::Path;
use treetime_utils::io::json::to_json_value;

pub fn load_config_document<T>(source: &ConfigSource, text: &str, base: Option<&Path>) -> Result<Value, Report>
where
  T: Serialize + Default + JsonSchema,
{
  let mut file_value = parse_config_document(source, text)?;
  if let Value::Object(map) = &mut file_value {
    map.remove(SCHEMA_KEY);
  }
  if let Some(base) = base {
    resolve_config_paths(&mut file_value, command_schema::<T>().as_value(), base)?;
  }
  let mut merged = to_json_value(&T::default())?;
  merge_value(&mut merged, &file_value);
  Ok(merged)
}

pub fn check_command_config(source: &ConfigSource, value: &Value, schema: &Schema) -> Result<(), Report> {
  let diags = schema_diagnostics(value, schema, false);
  render_and_bail(source, "invalid configuration", diags)
}

pub fn merge_value(base: &mut Value, overlay: &Value) {
  match (base, overlay) {
    (Value::Object(base), Value::Object(overlay)) => {
      for (key, value) in overlay {
        merge_value(base.entry(key.clone()).or_insert(Value::Null), value);
      }
    },
    (base, overlay) => *base = overlay.clone(),
  }
}
