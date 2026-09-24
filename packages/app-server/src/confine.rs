use app_commands::command::AppCommand;
use app_commands::commands::ancestral::aa_node_data::translation_input_paths;
use app_commands::config::properties::{PathRole, leaf_properties};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde_json::{Map, Value};
use std::path::{Path, PathBuf};
use treetime_utils::{make_error, make_report};

const OUTPUT_ALL_KEY: &str = "output_all";

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

  pub(crate) fn confine(&self, command: AppCommand, config: &mut Value, output_dir: &Path) -> Result<(), Report> {
    let Value::Object(settings) = config else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    let leaves = leaf_properties(command.config_schema().as_value())?;

    for leaf in &leaves {
      match leaf.path_role {
        Some(PathRole::Output) => remove_setting(settings, &leaf.key_path),
        Some(PathRole::Input) => {
          if let Some(value) = setting_mut(settings, &leaf.key_path) {
            self.confine_value(&leaf.key_path.join("."), value)?;
          }
        },
        Some(PathRole::InputTemplate) | None => {},
      }
    }

    for leaf in leaves
      .iter()
      .filter(|leaf| leaf.path_role == Some(PathRole::InputTemplate))
    {
      self.check_template(&leaf.key_path.join("."), settings, &leaf.key_path)?;
    }

    settings.insert(
      OUTPUT_ALL_KEY.to_owned(),
      Value::String(output_dir.to_string_lossy().into_owned()),
    );
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

  fn confine_path(&self, setting: &str, path: &Path) -> Result<PathBuf, Report> {
    let joined = self.data_dir.join(path);
    let resolved = joined.canonicalize().map_err(|err| {
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

  fn check_template(&self, setting: &str, settings: &Map<String, Value>, key_path: &[String]) -> Result<(), Report> {
    let Some(Value::String(template)) = setting_ref(settings, key_path) else {
      return Ok(());
    };
    let cdses: Vec<String> = settings
      .get("cdses")
      .and_then(Value::as_array)
      .map(|names| names.iter().filter_map(Value::as_str).map(str::to_owned).collect())
      .unwrap_or_default();
    let annotation = settings.get("annotation").and_then(Value::as_str).map(Path::new);
    let template = self.data_dir.join(template).to_string_lossy().into_owned();
    let paths = translation_input_paths(&template, &cdses, annotation)
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

fn setting_mut<'a>(settings: &'a mut Map<String, Value>, key_path: &[String]) -> Option<&'a mut Value> {
  let (first, rest) = key_path.split_first()?;
  rest
    .iter()
    .try_fold(settings.get_mut(first)?, |value, key| value.get_mut(key))
}

fn setting_ref<'a>(settings: &'a Map<String, Value>, key_path: &[String]) -> Option<&'a Value> {
  let (first, rest) = key_path.split_first()?;
  rest.iter().try_fold(settings.get(first)?, |value, key| value.get(key))
}

fn remove_setting(settings: &mut Map<String, Value>, key_path: &[String]) {
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
