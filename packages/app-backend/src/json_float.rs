use deser::adapters::{DeserializeAs, TryFromInto};
use deser::de::SinkHandle;
use deser::ser::Chunk;
use deser::{Atom, Deserialize, Error, Serialize, State};
use schemars::{JsonSchema, Schema, SchemaGenerator, json_schema};
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

#[derive(Deserialize)]
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

impl Serialize for JsonFloat {
  fn serialize(&self, _state: &mut State) -> Result<Chunk<'_>, Error> {
    let value = self.0;
    let atom = if value.is_finite() {
      Atom::F64(value)
    } else if value.is_nan() {
      Atom::Str(NOT_A_NUMBER.into())
    } else if value.is_sign_positive() {
      Atom::Str(POSITIVE_INFINITY.into())
    } else {
      Atom::Str(NEGATIVE_INFINITY.into())
    };
    Ok(Chunk::Atom(atom))
  }
}

impl<'de> Deserialize<'de> for JsonFloat {
  fn deserialize_into<'out>(out: &'out mut Option<Self>, state: &mut State) -> SinkHandle<'out, 'de> {
    <TryFromInto<JsonFloatRepr> as DeserializeAs<'de, Self>>::deserialize_into_as(out, state)
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
