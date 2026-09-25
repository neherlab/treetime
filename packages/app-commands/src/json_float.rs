use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde::de::Error as _;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Cow;

const POSITIVE_INFINITY: &str = "inf";
const NEGATIVE_INFINITY: &str = "-inf";
const NOT_A_NUMBER: &str = "nan";

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct JsonFloat(pub f64);

impl Serialize for JsonFloat {
  fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
    let value = self.0;
    if value.is_finite() {
      serializer.serialize_f64(value)
    } else if value.is_nan() {
      serializer.serialize_str(NOT_A_NUMBER)
    } else if value.is_sign_positive() {
      serializer.serialize_str(POSITIVE_INFINITY)
    } else {
      serializer.serialize_str(NEGATIVE_INFINITY)
    }
  }
}

impl<'de> Deserialize<'de> for JsonFloat {
  fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
    match JsonFloatRepr::deserialize(deserializer)? {
      JsonFloatRepr::Number(value) => Ok(Self(value)),
      JsonFloatRepr::Text(text) => match text.as_str() {
        POSITIVE_INFINITY => Ok(Self(f64::INFINITY)),
        NEGATIVE_INFINITY => Ok(Self(f64::NEG_INFINITY)),
        NOT_A_NUMBER => Ok(Self(f64::NAN)),
        other => Err(D::Error::custom(format!(
          "expected a number, \"{POSITIVE_INFINITY}\", \"{NEGATIVE_INFINITY}\" or \"{NOT_A_NUMBER}\", found \"{other}\""
        ))),
      },
    }
  }
}

#[derive(Deserialize)]
#[serde(untagged)]
enum JsonFloatRepr {
  Number(f64),
  Text(String),
}

impl JsonSchema for JsonFloat {
  fn schema_name() -> Cow<'static, str> {
    Cow::Borrowed("JsonFloat")
  }

  fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
    json_schema!({
      "description": "A number; infinities and NaN are the strings \"inf\", \"-inf\" and \"nan\".",
      "anyOf": [
        { "type": "number" },
        { "type": "string", "enum": [POSITIVE_INFINITY, NEGATIVE_INFINITY, NOT_A_NUMBER] }
      ]
    })
  }
}
