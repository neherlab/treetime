use schemars::generate::SchemaSettings;
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde_json::json;

pub const SCHEMA_KEY: &str = "$schema";

pub fn command_schema<T: JsonSchema>() -> Schema {
  let mut schema = draft2020_generator().into_root_schema_for::<T>();
  allow_schema_ref(&mut schema);
  schema
}

pub fn draft2020_generator() -> SchemaGenerator {
  draft2020_settings().into_generator()
}

pub fn draft2020_settings() -> SchemaSettings {
  SchemaSettings::draft2020_12()
}

fn allow_schema_ref(schema: &mut Schema) {
  let object = schema.ensure_object();
  let properties = object.entry("properties").or_insert_with(|| json!({}));
  if let Some(properties) = properties.as_object_mut() {
    properties.insert(
      SCHEMA_KEY.to_owned(),
      json!({
        "type": "string",
        "description": "Path or URL of the JSON schema for this config; used by editors and ignored by the loader."
      }),
    );
  }
}
