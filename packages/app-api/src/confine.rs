use app_backend::runs::inputs::template_input_paths;
use app_commands::command::AppCommand;
use app_commands::config::properties::{PathRole, setting_properties};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde_json::{Map, Value};
use std::path::{Path, PathBuf};
use treetime_utils::io::fs::absolute_path;
use treetime_utils::{make_error, make_report};

#[derive(Clone, Debug)]
pub(crate) struct PathPolicy {
  input_roots: Vec<PathBuf>,
}

impl PathPolicy {
  pub(crate) fn new(input_roots: &[PathBuf]) -> Result<Self, Report> {
    let input_roots = input_roots
      .iter()
      .map(|root| canonical_dir(root).wrap_err("When resolving an input directory"))
      .try_collect()?;
    Ok(Self { input_roots })
  }

  pub(crate) fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report> {
    let Value::Object(settings) = config else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    let properties = setting_properties(command.config_schema().as_value())?;

    for property in &properties {
      match property.path_role {
        Some(PathRole::Input) => {
          if let Some(value) = settings.get_mut(&property.key) {
            self.confine_value(&property.key, value)?;
          }
        },
        Some(PathRole::InputTemplate | PathRole::Output) | None => {},
      }
    }

    for property in properties
      .iter()
      .filter(|property| property.path_role == Some(PathRole::InputTemplate))
    {
      self.confine_template(&property.key, settings)?;
    }

    Ok(())
  }

  fn confine_value(&self, setting: &str, value: &mut Value) -> Result<(), Report> {
    match value {
      Value::String(path) => {
        *path = self
          .confine_path(setting, Path::new(path))?
          .to_string_lossy()
          .into_owned();
      },
      Value::Array(items) => {
        for item in items {
          if let Value::String(path) = item {
            *path = self
              .confine_path(setting, Path::new(path))?
              .to_string_lossy()
              .into_owned();
          }
        }
      },
      _ => {},
    }
    Ok(())
  }

  pub(crate) fn confine_path(&self, setting: &str, path: &Path) -> Result<PathBuf, Report> {
    let resolved = path.canonicalize().map_err(|err| {
      make_report!(
        "input `{}` of setting `{setting}` cannot be read: {err}",
        path.display()
      )
    })?;
    if self.input_roots.iter().any(|root| resolved.starts_with(root)) {
      Ok(resolved)
    } else {
      make_error!(
        "input `{}` of setting `{setting}` is outside the directories the server reads inputs from",
        path.display()
      )
    }
  }

  fn confine_template(&self, setting: &str, settings: &mut Map<String, Value>) -> Result<(), Report> {
    let Some(Value::String(template)) = settings.get_mut(setting) else {
      return Ok(());
    };
    *template = absolute_path(&*template)?.to_string_lossy().into_owned();

    let paths = template_input_paths(settings, setting)
      .wrap_err_with(|| format!("When listing the inputs of setting `{setting}`"))?;
    let rejected = paths
      .iter()
      .filter_map(|path| self.confine_path(setting, path).err())
      .map(|err| err.to_string())
      .collect_vec();
    if rejected.is_empty() {
      Ok(())
    } else {
      make_error!("{}", rejected.join("; "))
    }
  }
}

fn canonical_dir(dir: &Path) -> Result<PathBuf, Report> {
  let resolved = dir
    .canonicalize()
    .wrap_err_with(|| format!("When resolving directory '{}'", dir.display()))?;
  if resolved.is_dir() {
    Ok(resolved)
  } else {
    make_error!("'{}' is not a directory", dir.display())
  }
}
