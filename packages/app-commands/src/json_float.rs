use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
use serde::de::Error as _;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::borrow::Cow;
use treetime_primitives::LogLh;

const POSITIVE_INFINITY: &str = "inf";
const NEGATIVE_INFINITY: &str = "-inf";
const NOT_A_NUMBER: &str = "nan";

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct JsonFloat(pub f64);

impl From<LogLh> for JsonFloat {
  fn from(value: LogLh) -> Self {
    Self(value.value())
  }
}

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

#[derive(Deserialize, deser::Deserialize)]
#[serde(untagged)]
#[deser(untagged)]
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

impl deser::Serialize for JsonFloat {
  fn serialize(&self, _state: &mut deser::State) -> Result<deser::ser::Chunk<'_>, deser::Error> {
    let value = self.0;
    let atom = if value.is_finite() {
      deser::Atom::F64(value)
    } else if value.is_nan() {
      deser::Atom::Str(NOT_A_NUMBER.into())
    } else if value.is_sign_positive() {
      deser::Atom::Str(POSITIVE_INFINITY.into())
    } else {
      deser::Atom::Str(NEGATIVE_INFINITY.into())
    };
    Ok(deser::ser::Chunk::Atom(atom))
  }
}

impl<'de> deser::Deserialize<'de> for JsonFloat {
  fn deserialize_into<'out>(out: &'out mut Option<Self>, state: &mut deser::State) -> deser::de::SinkHandle<'out, 'de> {
    <deser::adapters::TryFromInto<JsonFloatRepr> as deser::adapters::DeserializeAs<'de, Self>>::deserialize_into_as(
      out, state,
    )
  }
}

impl TryFrom<JsonFloatRepr> for JsonFloat {
  type Error = String;

  fn try_from(repr: JsonFloatRepr) -> Result<Self, String> {
    match repr {
      JsonFloatRepr::Number(value) => Ok(Self(value)),
      JsonFloatRepr::Text(text) => match text.as_str() {
        POSITIVE_INFINITY => Ok(Self(f64::INFINITY)),
        NEGATIVE_INFINITY => Ok(Self(f64::NEG_INFINITY)),
        NOT_A_NUMBER => Ok(Self(f64::NAN)),
        other => Err(format!(
          "expected a number, \"{POSITIVE_INFINITY}\", \"{NEGATIVE_INFINITY}\" or \"{NOT_A_NUMBER}\", found \"{other}\""
        )),
      },
    }
  }
}
