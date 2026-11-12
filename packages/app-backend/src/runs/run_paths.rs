use crate::runs::errors::invalid;
use crate::runs::input_store::{InputStore, StoredInput};
use crate::runs::record::RunInput;
use app_commands::command::AppCommand;
use app_commands::config::properties::setting_properties;
use eyre::Report;
use itertools::Itertools;
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

pub const INPUTS_DIR: &str = "inputs";

pub const OUT_DIR: &str = "out";

pub fn resolve_settings(
  command: AppCommand,
  settings: &Map<String, Value>,
  inputs: &[RunInput],
  run_dir: &Path,
  store: &InputStore,
) -> Result<Map<String, Value>, Report> {
  let mut settings = settings.clone();
  map_path_values(command, &mut settings, |path| {
    resolve_path(path, inputs, run_dir, store)
  })?;
  Ok(settings)
}

pub fn resolve_path(path: &Path, inputs: &[RunInput], run_dir: &Path, store: &InputStore) -> PathBuf {
  if path.is_absolute() {
    return path.to_path_buf();
  }
  inputs
    .iter()
    .find(|input| input.path == path)
    .and_then(|input| input.stored())
    .map_or_else(|| run_dir.join(path), |key| store.path(&key))
}

pub struct PortablePaths {
  run_dir: Option<PathBuf>,
  stored: BTreeMap<PathBuf, StoredPath>,
}

impl PortablePaths {
  pub fn of_settings(command: AppCommand, settings: &Map<String, Value>, store: &InputStore) -> Result<Self, Report> {
    let mut found: Vec<(PathBuf, StoredInput)> = vec![];
    let mut settings = settings.clone();
    map_path_values(command, &mut settings, |path| {
      if let Some(key) = store.key_of(path)
        && !found.iter().any(|(known, _)| known == path)
      {
        found.push((path.to_path_buf(), key));
      }
      path.to_path_buf()
    })?;
    let by_name = found.iter().into_group_map_by(|(_, key)| key.name.clone());
    let stored = found
      .iter()
      .map(|(path, key)| {
        let clashing = &by_name[&key.name];
        let portable = if clashing.len() == 1 {
          Path::new(INPUTS_DIR).join(&key.name)
        } else {
          let number = clashing
            .iter()
            .position(|(other, _)| other == path)
            .map_or(0, |index| index + 1);
          Path::new(INPUTS_DIR).join(number.to_string()).join(&key.name)
        };
        let size = store.size(key).map_err(|report| {
          invalid(format!(
            "input `{}` is not a file uploaded to TreeTime: {report:#}",
            path.display()
          ))
        })?;
        Ok((
          path.clone(),
          StoredPath {
            portable,
            key: key.clone(),
            size,
          },
        ))
      })
      .collect::<Result<BTreeMap<_, _>, Report>>()?;
    Ok(Self { run_dir: None, stored })
  }

  #[must_use]
  pub fn in_run(self, run_dir: &Path) -> Self {
    Self {
      run_dir: Some(run_dir.to_path_buf()),
      ..self
    }
  }

  pub fn portable(&self, path: &Path) -> PathBuf {
    if let Some(stored) = self.stored.get(path) {
      return stored.portable.clone();
    }
    self
      .run_dir
      .as_ref()
      .and_then(|run_dir| path.strip_prefix(run_dir).ok())
      .map_or_else(|| path.to_path_buf(), Path::to_path_buf)
  }

  pub fn apply(&self, command: AppCommand, settings: &mut Map<String, Value>) -> Result<(), Report> {
    map_path_values(command, settings, |path| self.portable(path))
  }

  pub fn stored_inputs(&self, command: AppCommand, settings: &Map<String, Value>) -> Result<Vec<RunInput>, Report> {
    let mut inputs = vec![];
    for property in setting_properties(command.config_schema().as_value())? {
      if property.path_role.is_none() {
        continue;
      }
      let Some(value) = settings.get(&property.key) else {
        continue;
      };
      for path in path_strings(value) {
        if let Some(stored) = self.stored.get(Path::new(path)) {
          inputs.push(RunInput {
            setting: property.key.clone(),
            path: stored.portable.clone(),
            size: stored.size,
            sha256: stored.key.sha256.clone(),
          });
        }
      }
    }
    Ok(inputs)
  }
}

struct StoredPath {
  portable: PathBuf,
  key: StoredInput,
  size: usize,
}

pub fn map_path_values(
  command: AppCommand,
  settings: &mut Map<String, Value>,
  mut map: impl FnMut(&Path) -> PathBuf,
) -> Result<(), Report> {
  for property in setting_properties(command.config_schema().as_value())? {
    if property.path_role.is_none() {
      continue;
    }
    if let Some(value) = settings.get_mut(&property.key) {
      map_strings(value, &mut map);
    }
  }
  Ok(())
}

fn map_strings(value: &mut Value, map: &mut impl FnMut(&Path) -> PathBuf) {
  match value {
    Value::String(path) if !path.is_empty() => {
      *path = map(Path::new(path.as_str())).to_string_lossy().into_owned();
    },
    Value::Array(items) => {
      for item in items {
        map_strings(item, map);
      }
    },
    _ => {},
  }
}

pub fn path_strings(value: &Value) -> Vec<&str> {
  match value {
    Value::String(path) if !path.is_empty() => vec![path.as_str()],
    Value::Array(items) => items.iter().flat_map(path_strings).collect(),
    _ => vec![],
  }
}
