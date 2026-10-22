use crate::config::resolve_paths::resolve_config_paths;
use crate::config::schema::{SCHEMA_KEY, command_schema};
use crate::config::schema_check::schema_diagnostics;
use crate::config::source::{
  ConfigSource, RawDiagnostic, escape_pointer, invalid_config_report, parse_config_document, render_and_bail,
};
use derive_more::{Display, Error};
use deser::Serialize;
use deser::adapters::As;
use deser::de::DeserializeOwned;
use deser_path::{PathLayer, PathSegment};
use deser_serde::Serde;
use deser_value::Deserializer;
use eyre::{Report, WrapErr};
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

pub fn read_settings<T: DeserializeOwned>(settings: &Value) -> Result<T, Report> {
  let value = deser_value::to_value(&As::<_, Serde>::new(settings)).wrap_err("When converting a JSON value")?;
  Deserializer::new(&value)
    .deserialize_with(|driver| driver.push_layer(PathLayer::new()))
    .map_err(|err| match err.attachment::<deser_path::Path>() {
      Some(path) => Report::new(InvalidSetting {
        pointer: setting_pointer(path.segments()),
        message: err.message().to_owned(),
      }),
      None => Report::new(err).wrap_err("When converting a JSON value"),
    })
}

pub fn invalid_setting_report(source: &ConfigSource, report: Report) -> Report {
  let Some(invalid) = report.chain().find_map(|error| error.downcast_ref::<InvalidSetting>()) else {
    return report;
  };
  invalid_config_report(
    source,
    "invalid configuration",
    vec![
      RawDiagnostic::builder("config::value", invalid.message.clone())
        .at(invalid.pointer.clone())
        .build(),
    ],
  )
}

#[derive(Debug, Display, Error)]
#[display("{message}")]
pub struct InvalidSetting {
  pub pointer: String,
  pub message: String,
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

fn setting_pointer(segments: &[PathSegment]) -> String {
  segments
    .iter()
    .map(|segment| match segment {
      PathSegment::Key(key) => format!("/{}", escape_pointer(key)),
      PathSegment::Index(index) => format!("/{index}"),
      PathSegment::Unknown => String::new(),
    })
    .collect()
}
