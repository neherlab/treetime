use crate::config::schema::SCHEMA_KEY;
use crate::config::source::escape_pointer;
use eyre::Report;
use serde_json::Value;
use std::iter;
use strum_macros::EnumString;
use treetime_utils::{make_error, make_report};

pub const PATH_ROLE_KEY: &str = "x-path";

pub const CLI_FLAG_KEY: &str = "x-cli-flag";

pub const DEFS_PREFIX: &str = "#/$defs/";

pub fn setting_properties(schema: &Value) -> Result<Vec<SettingProperty>, Report> {
  let Some(properties) = schema.get("properties").and_then(Value::as_object) else {
    return Ok(vec![]);
  };
  properties
    .iter()
    .filter(|(key, _)| *key != SCHEMA_KEY)
    .map(|(key, property)| {
      let schema_pointer = format!("/properties/{}", escape_pointer(key));
      if is_setting_group(schema, property) {
        return make_error!(
          "config key `{key}` holds a group of settings: every setting of a command must be a top-level key of its config"
        );
      }
      let path_role = property
        .get(PATH_ROLE_KEY)
        .and_then(Value::as_str)
        .map(|role| {
          role
            .parse()
            .map_err(|err| make_report!("unknown `{PATH_ROLE_KEY}` value `{role}` at `{schema_pointer}`: {err}"))
        })
        .transpose()?;
      Ok(SettingProperty {
        key: key.clone(),
        schema_pointer,
        path_role,
      })
    })
    .collect()
}

pub struct SettingProperty {
  pub key: String,
  pub schema_pointer: String,
  pub path_role: Option<PathRole>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, EnumString)]
#[strum(serialize_all = "kebab-case")]
pub enum PathRole {
  Input,
  InputTemplate,
  Output,
}

fn is_setting_group(schema: &Value, property: &Value) -> bool {
  iter::once(property)
    .chain(
      ["anyOf", "oneOf", "allOf"]
        .iter()
        .filter_map(|key| property.get(*key))
        .filter_map(Value::as_array)
        .flatten(),
    )
    .any(|candidate| {
      candidate.get("properties").is_some()
        || candidate
          .get("$ref")
          .and_then(Value::as_str)
          .and_then(|reference| reference.strip_prefix(DEFS_PREFIX))
          .and_then(|def| schema.pointer(&def_pointer(def)))
          .is_some_and(|def| def.get("properties").is_some())
    })
}

pub fn def_pointer(name: &str) -> String {
  format!("/$defs/{}", escape_pointer(name))
}
