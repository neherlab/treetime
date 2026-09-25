use app_commands::command::AppCommand;
use app_commands::config::properties::{PathRole, leaf_properties};
use app_commands::config::settings::setting_mut;
use app_commands::runs::inputs::template_input_paths;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde_json::{Map, Value};
use std::path::{Path, PathBuf};
use treetime_utils::io::fs::absolute_path;
use treetime_utils::{make_error, make_report};

#[derive(Clone, Debug)]
pub(crate) struct PathPolicy {
  data_dir: PathBuf,
  input_roots: Vec<PathBuf>,
}

impl PathPolicy {
  pub(crate) fn new(data_dir: &Path, extra_input_roots: &[PathBuf]) -> Result<Self, Report> {
    let data_dir = canonical_dir(data_dir).wrap_err("When resolving the data directory")?;
    let extra_input_roots: Vec<PathBuf> = extra_input_roots
      .iter()
      .map(|root| canonical_dir(root).wrap_err("When resolving an input directory"))
      .try_collect()?;
    Ok(Self {
      input_roots: [vec![data_dir.clone()], extra_input_roots].concat(),
      data_dir,
    })
  }

  pub(crate) fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report> {
    let Value::Object(settings) = config else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    let leaves = leaf_properties(command.config_schema().as_value())?;

    for leaf in &leaves {
      match leaf.path_role {
        Some(PathRole::Input) => {
          if let Some(value) = setting_mut(settings, &leaf.key_path) {
            self.confine_value(&leaf.key_path.join("."), value)?;
          }
        },
        Some(PathRole::InputTemplate | PathRole::Output) | None => {},
      }
    }

    for leaf in leaves
      .iter()
      .filter(|leaf| leaf.path_role == Some(PathRole::InputTemplate))
    {
      self.confine_template(&leaf.key_path.join("."), settings, &leaf.key_path)?;
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
    let in_data_dir = self.data_dir.join(path);
    let candidate = if path.is_relative() && !in_data_dir.exists() {
      path.to_path_buf()
    } else {
      in_data_dir
    };
    let resolved = candidate.canonicalize().map_err(|err| {
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

  fn confine_template(
    &self,
    setting: &str,
    settings: &mut Map<String, Value>,
    key_path: &[String],
  ) -> Result<(), Report> {
    let Some(Value::String(template)) = setting_mut(settings, key_path) else {
      return Ok(());
    };
    let in_data_dir = self.data_dir.join(&*template);
    *template = if Path::new(template).is_relative() && !in_data_dir.parent().is_some_and(Path::exists) {
      absolute_path(&*template)?
    } else {
      in_data_dir
    }
    .to_string_lossy()
    .into_owned();

    let paths = template_input_paths(settings, key_path)
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
