use crate::config::schema::SCHEMA_KEY;
use crate::config::source::escape_pointer;
use eyre::Report;
use serde_json::Value;
use std::collections::BTreeSet;
use std::{iter, slice};
use strum_macros::EnumString;
use treetime_utils::{make_error, make_report};

pub const PATH_ROLE_KEY: &str = "x-path";

pub const CLI_FLAG_KEY: &str = "x-cli-flag";

pub const CLI_NUM_ARGS_KEY: &str = "x-cli-num-args";

pub const CLI_VALUE_DELIMITER_KEY: &str = "x-cli-value-delimiter";

pub const CLI_VALUES_KEY: &str = "x-cli-values";

pub fn leaf_properties(schema: &Value) -> Result<Vec<LeafProperty>, Report> {
  let mut leaves = Vec::new();
  let mut visited_defs = BTreeSet::new();
  collect_leaves(schema, "", &[], &mut visited_defs, &mut leaves)?;
  Ok(leaves)
}

pub struct LeafProperty {
  pub key_path: Vec<String>,
  pub schema_pointer: String,
  pub path_role: Option<PathRole>,
}

impl LeafProperty {
  pub fn key(&self) -> &str {
    self.key_path.last().map_or("", String::as_str)
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, EnumString)]
#[strum(serialize_all = "kebab-case")]
pub enum PathRole {
  Input,
  InputTemplate,
  Output,
}

fn collect_leaves(
  schema: &Value,
  object_pointer: &str,
  key_path: &[String],
  visited_defs: &mut BTreeSet<String>,
  leaves: &mut Vec<LeafProperty>,
) -> Result<(), Report> {
  let Some(properties) = schema
    .pointer(object_pointer)
    .and_then(|object| object.get("properties"))
    .and_then(Value::as_object)
  else {
    return Ok(());
  };
  for (key, property) in properties {
    if key_path.is_empty() && key == SCHEMA_KEY {
      continue;
    }
    let key_path = [key_path, slice::from_ref(key)].concat();
    let schema_pointer = format!("{object_pointer}/properties/{}", escape_pointer(key));
    if let Some(def) = nested_object_def(schema, property) {
      if !visited_defs.insert(def.clone()) {
        return make_error!(
          "schema definition `{def}` holds the settings of more than one config key, so its settings cannot carry per-key annotations"
        );
      }
      collect_leaves(
        schema,
        &format!("/$defs/{}", escape_pointer(&def)),
        &key_path,
        visited_defs,
        leaves,
      )?;
    } else {
      let path_role = match property.get(PATH_ROLE_KEY).and_then(Value::as_str) {
        Some(role) => Some(
          role
            .parse()
            .map_err(|err| make_report!("unknown `{PATH_ROLE_KEY}` value `{role}` at `{schema_pointer}`: {err}"))?,
        ),
        None => None,
      };
      leaves.push(LeafProperty {
        key_path,
        schema_pointer,
        path_role,
      });
    }
  }
  Ok(())
}

fn nested_object_def(schema: &Value, property: &Value) -> Option<String> {
  let refs = iter::once(property).chain(
    ["anyOf", "oneOf", "allOf"]
      .iter()
      .filter_map(|key| property.get(*key))
      .filter_map(Value::as_array)
      .flatten(),
  );
  refs
    .filter_map(|candidate| candidate.get("$ref").and_then(Value::as_str))
    .filter_map(|reference| reference.strip_prefix("#/$defs/"))
    .find(|def| {
      schema
        .pointer(&format!("/$defs/{}", escape_pointer(def)))
        .is_some_and(|def| def.get("properties").is_some())
    })
    .map(str::to_owned)
}
