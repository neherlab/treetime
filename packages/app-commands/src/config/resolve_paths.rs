use crate::config::properties::{PathRole, leaf_properties};
use crate::config::settings::setting_mut;
use eyre::Report;
use serde_json::Value;
use std::path::Path;
use treetime_utils::make_report;

const STDIO_PATH: &str = "-";

pub fn resolve_config_paths(config: &mut Value, schema: &Value, base: &Path) -> Result<(), Report> {
  resolve_config_paths_where(config, schema, base, |_| true)
}

pub fn resolve_config_paths_where(
  config: &mut Value,
  schema: &Value,
  base: &Path,
  include: impl Fn(&[String]) -> bool,
) -> Result<(), Report> {
  let Value::Object(settings) = config else {
    return Ok(());
  };
  for leaf in leaf_properties(schema)? {
    if !matches!(
      leaf.path_role,
      Some(PathRole::Input | PathRole::InputTemplate | PathRole::Output)
    ) || !include(&leaf.key_path)
    {
      continue;
    }
    if let Some(value) = setting_mut(settings, &leaf.key_path) {
      resolve_value(value, base)?;
    }
  }
  Ok(())
}

fn resolve_value(value: &mut Value, base: &Path) -> Result<(), Report> {
  match value {
    Value::String(path) => resolve_string(path, base),
    Value::Array(paths) => paths.iter_mut().try_for_each(|path| resolve_value(path, base)),
    _ => Ok(()),
  }
}

fn resolve_string(path: &mut String, base: &Path) -> Result<(), Report> {
  if path.is_empty() || path == STDIO_PATH || Path::new(path.as_str()).is_absolute() {
    return Ok(());
  }
  let joined = base.join(path.as_str());
  *path = joined
    .to_str()
    .ok_or_else(|| make_report!("the path '{}' is not valid UTF-8", joined.display()))?
    .to_owned();
  Ok(())
}
