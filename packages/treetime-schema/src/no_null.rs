use schemars::Schema;
use schemars::transform::{Transform, transform_subschemas};
use serde_json::{Map, Value};

pub const UNSET_KEY: &str = "x-unset";

const NULL_TYPE: &str = "null";

#[derive(Clone, Copy, Debug, Default)]
pub struct NoNull;

impl Transform for NoNull {
  fn transform(&mut self, schema: &mut Schema) {
    transform_subschemas(self, schema);
    let Some(object) = schema.as_object_mut() else {
      return;
    };
    let removed = [
      remove_null_type(object),
      remove_null_branch(object),
      remove_null_enum_value(object),
      remove_null_default(object),
    ]
    .contains(&true);
    if removed {
      object.insert(UNSET_KEY.to_owned(), Value::Bool(true));
    }
  }
}

fn remove_null_type(object: &mut Map<String, Value>) -> bool {
  let Some(Value::Array(types)) = object.get_mut("type") else {
    return false;
  };
  let before = types.len();
  types.retain(|ty| ty != NULL_TYPE);
  if types.len() == before {
    return false;
  }
  if let [single] = types.as_slice() {
    let single = single.clone();
    object.insert("type".to_owned(), single);
  }
  true
}

fn remove_null_branch(object: &mut Map<String, Value>) -> bool {
  let Some(Value::Array(branches)) = object.get_mut("anyOf") else {
    return false;
  };
  let before = branches.len();
  branches.retain(|branch| !is_null_schema(branch));
  if branches.len() == before {
    return false;
  }
  if let [Value::Object(single)] = branches.as_slice() {
    let single = single.clone();
    object.remove("anyOf");
    for (key, value) in single {
      object.entry(key).or_insert(value);
    }
  }
  true
}

fn remove_null_enum_value(object: &mut Map<String, Value>) -> bool {
  let Some(Value::Array(values)) = object.get_mut("enum") else {
    return false;
  };
  let before = values.len();
  values.retain(|value| !value.is_null());
  if values.len() == before {
    return false;
  }
  if let [single] = values.as_slice() {
    let single = single.clone();
    object.remove("enum");
    object.insert("const".to_owned(), single);
  }
  true
}

fn remove_null_default(object: &mut Map<String, Value>) -> bool {
  if object.get("default").is_some_and(Value::is_null) {
    object.remove("default");
    true
  } else {
    false
  }
}

fn is_null_schema(schema: &Value) -> bool {
  schema
    .as_object()
    .is_some_and(|object| object.len() == 1 && object.get("type").is_some_and(|ty| ty == NULL_TYPE))
}
