use serde_json::{Map, Value};

pub fn setting_ref<'a>(settings: &'a Map<String, Value>, key_path: &[String]) -> Option<&'a Value> {
  let (first, rest) = key_path.split_first()?;
  rest.iter().try_fold(settings.get(first)?, |value, key| value.get(key))
}

pub fn setting_mut<'a>(settings: &'a mut Map<String, Value>, key_path: &[String]) -> Option<&'a mut Value> {
  let (first, rest) = key_path.split_first()?;
  rest
    .iter()
    .try_fold(settings.get_mut(first)?, |value, key| value.get_mut(key))
}

pub fn remove_setting(settings: &mut Map<String, Value>, key_path: &[String]) {
  let Some((last, parents)) = key_path.split_last() else {
    return;
  };
  let parent = if parents.is_empty() {
    Some(settings)
  } else {
    setting_mut(settings, parents).and_then(Value::as_object_mut)
  };
  if let Some(parent) = parent {
    parent.remove(last);
  }
}

pub fn has_path(value: &Value) -> bool {
  match value {
    Value::String(path) => !path.is_empty(),
    Value::Array(paths) => paths.iter().any(has_path),
    _ => false,
  }
}
