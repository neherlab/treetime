use crate::config::schema::SCHEMA_KEY;
use crate::config::schema_check::schema_diagnostics;
use crate::config::source::{ConfigSource, parse_config_document, render_and_bail};
use eyre::Report;
use schemars::Schema;
use serde::Serialize;
use serde_json::Value;

pub fn load_config_document<T>(source: &ConfigSource, text: &str) -> Result<Value, Report>
where
  T: Serialize + Default,
{
  let mut file_value = parse_config_document(source, text)?;
  if let Value::Object(map) = &mut file_value {
    map.remove(SCHEMA_KEY);
  }
  let mut merged = serde_json::to_value(T::default())?;
  merge_value(&mut merged, &file_value);
  Ok(merged)
}

pub fn check_command_config(source: &ConfigSource, value: &Value, schema: &Schema) -> Result<(), Report> {
  let diags = schema_diagnostics(value, schema, false);
  render_and_bail(source, "invalid configuration", diags)
}

fn merge_value(base: &mut Value, overlay: &Value) {
  match (base, overlay) {
    (Value::Object(base), Value::Object(overlay)) => {
      for (key, value) in overlay {
        merge_value(base.entry(key.clone()).or_insert(Value::Null), value);
      }
    },
    (base, overlay) => *base = overlay.clone(),
  }
}
