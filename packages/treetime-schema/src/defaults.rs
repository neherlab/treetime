use deser::Serialize;
use schemars::Schema;
use serde_json::{Map, Value};
use std::mem::take;
use treetime_utils::io::json::to_json_value;

const EXTENSION_PREFIX: &str = "x-";

#[expect(
  clippy::expect_used,
  reason = "a schema transform cannot return errors, and the default of a configuration type always serializes"
)]
pub fn schema_defaults<T: Default + Serialize>(schema: &mut Schema) {
  let Value::Object(defaults) = to_json_value(&T::default()).expect("the default value serializes to JSON") else {
    return;
  };
  let Some(Value::Object(properties)) = schema.get_mut("properties") else {
    return;
  };
  for (key, default) in defaults {
    if let Some(Value::Object(property)) = properties.get_mut(&key)
      && !property.contains_key("default")
    {
      insert_before_extensions(property, default);
    }
  }
}

#[expect(
  clippy::expect_used,
  reason = "a schema transform cannot return errors, and the default of a configuration type always serializes"
)]
pub fn field_default<T: Default + Serialize>(schema: &mut Schema) {
  let default = to_json_value(&T::default()).expect("the default value serializes to JSON");
  if let Some(property) = schema.as_object_mut()
    && !property.contains_key("default")
  {
    insert_before_extensions(property, default);
  }
}

pub fn skip_serializing_optionals(schema: &mut Schema) {
  if let Some(object) = schema.as_object_mut() {
    unrequire_nullable(object);
  }
}

fn unrequire_nullable(object: &mut Map<String, Value>) {
  for key in ["oneOf", "anyOf"] {
    if let Some(Value::Array(branches)) = object.get_mut(key) {
      branches
        .iter_mut()
        .filter_map(Value::as_object_mut)
        .for_each(unrequire_nullable);
    }
  }
  let nullable = object
    .get("properties")
    .and_then(Value::as_object)
    .map(|properties| {
      properties
        .iter()
        .filter(|(_, property)| allows_null(property))
        .map(|(name, _)| Value::String(name.clone()))
        .collect::<Vec<_>>()
    })
    .unwrap_or_default();
  if let Some(Value::Array(required)) = object.get_mut("required") {
    required.retain(|name| !nullable.contains(name));
  }
}

fn allows_null(schema: &Value) -> bool {
  let is_null = |ty: &Value| ty == "null";
  match schema.get("type") {
    Some(Value::Array(types)) if types.iter().any(is_null) => return true,
    Some(ty) if is_null(ty) => return true,
    _ => {},
  }
  ["anyOf", "oneOf"]
    .iter()
    .filter_map(|key| schema.get(key).and_then(Value::as_array))
    .any(|branches| branches.iter().any(allows_null))
}

fn insert_before_extensions(property: &mut Map<String, Value>, default: Value) {
  let entries = take(property);
  let (head, tail): (Vec<_>, Vec<_>) = entries
    .into_iter()
    .partition(|(key, _)| !key.starts_with(EXTENSION_PREFIX));
  property.extend(head);
  property.insert("default".to_owned(), default);
  property.extend(tail);
}
