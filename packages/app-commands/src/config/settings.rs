use serde_json::Value;

pub fn has_path(value: &Value) -> bool {
  match value {
    Value::String(path) => !path.is_empty(),
    Value::Array(paths) => paths.iter().any(has_path),
    _ => false,
  }
}
