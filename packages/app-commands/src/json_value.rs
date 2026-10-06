use derive_more::{Deref, DerefMut, From, Into};
use deser_serde::Serde;
use schemars::generate::SchemaGenerator;
use schemars::{JsonSchema, Schema, json_schema};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::borrow::Cow;

#[derive(
  Clone,
  Debug,
  Default,
  PartialEq,
  Eq,
  Serialize,
  Deserialize,
  Deref,
  DerefMut,
  From,
  Into,
  deser::Serialize,
  deser::Deserialize,
)]
#[serde(transparent)]
pub struct JsonValue(#[deser(as = Serde)] pub Value);

impl JsonSchema for JsonValue {
  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("JsonValue")
  }

  fn json_schema(generator: &mut SchemaGenerator) -> Schema {
    let value = generator.subschema_for::<Self>();
    json_schema!({
      "description": "A JSON value other than `null`: a boolean, a number, a string, a list, or a mapping of them.",
      "anyOf": [
        { "type": "boolean" },
        { "type": "number" },
        { "type": "string" },
        { "type": "array", "items": value },
        { "type": "object", "additionalProperties": value },
      ],
    })
  }
}

#[derive(
  Clone,
  Debug,
  Default,
  PartialEq,
  Eq,
  Serialize,
  Deserialize,
  Deref,
  DerefMut,
  From,
  Into,
  deser::Serialize,
  deser::Deserialize,
)]
#[serde(transparent)]
pub struct SparseConfig(#[deser(as = Serde)] pub Map<String, Value>);

impl SparseConfig {
  pub fn into_value(self) -> Value {
    Value::Object(self.0)
  }
}

impl JsonSchema for SparseConfig {
  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("SparseConfig")
  }

  fn json_schema(generator: &mut SchemaGenerator) -> Schema {
    let value = generator.subschema_for::<JsonValue>();
    json_schema!({
      "description": "Settings of a command by key, as in a configuration file. Settings that are not listed take \
        their defaults.",
      "type": "object",
      "additionalProperties": value,
    })
  }
}
