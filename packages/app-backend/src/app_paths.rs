use crate::app_settings::settings::AppPathSettings;
use dirs::config_local_dir;
use eyre::Report;
use std::path::{Path, PathBuf};
use treetime_utils::env::env_var_optional;
use treetime_utils::io::fs::absolute_path;
use treetime_utils::make_report;

pub const APP_DIR_ENV: &str = "TREETIME_APP_DIR";

const APP_NAME: &str = "treetime";

pub fn app_root() -> Result<PathBuf, Report> {
  match env_path(APP_DIR_ENV)? {
    Some(root) => Ok(root),
    None => Ok(
      config_local_dir()
        .ok_or_else(|| make_report!("the configuration directory of the user could not be determined"))?
        .join(APP_NAME),
    ),
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AppPaths {
  pub root: PathBuf,
  pub profile: AppFolderPath,
  pub runs: AppFolderPath,
  pub logs: AppFolderPath,
  pub examples: AppFolderPath,
}

impl AppPaths {
  pub fn resolve(root: &Path, env: &AppFolderEnv, settings: &AppPathSettings) -> Self {
    let folder = |folder: AppFolder, from_env: Option<&PathBuf>, from_settings: Option<&PathBuf>| match from_env {
      Some(path) => AppFolderPath {
        path: path.clone(),
        fixed_by: Some(folder.env_var()),
      },
      None => AppFolderPath {
        path: from_settings.map_or_else(|| root.join(folder.dir_name()), |path| root.join(path)),
        fixed_by: None,
      },
    };
    Self {
      root: root.to_path_buf(),
      profile: folder(AppFolder::Profile, env.profile.as_ref(), settings.profile.as_ref()),
      runs: folder(AppFolder::Runs, env.runs.as_ref(), settings.runs.as_ref()),
      logs: folder(AppFolder::Logs, env.logs.as_ref(), settings.logs.as_ref()),
      examples: folder(AppFolder::Examples, env.examples.as_ref(), settings.examples.as_ref()),
    }
  }

  pub fn default_runs(&self) -> PathBuf {
    self.root.join(AppFolder::Runs.dir_name())
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AppFolderPath {
  pub path: PathBuf,
  pub fixed_by: Option<&'static str>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AppFolderEnv {
  pub profile: Option<PathBuf>,
  pub runs: Option<PathBuf>,
  pub logs: Option<PathBuf>,
  pub examples: Option<PathBuf>,
}

impl AppFolderEnv {
  pub fn from_env() -> Result<Self, Report> {
    Ok(Self {
      profile: env_path(AppFolder::Profile.env_var())?,
      runs: env_path(AppFolder::Runs.env_var())?,
      logs: env_path(AppFolder::Logs.env_var())?,
      examples: env_path(AppFolder::Examples.env_var())?,
    })
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AppFolder {
  Profile,
  Runs,
  Logs,
  Examples,
}

impl AppFolder {
  pub const fn dir_name(self) -> &'static str {
    match self {
      Self::Profile => "profile",
      Self::Runs => "runs",
      Self::Logs => "logs",
      Self::Examples => "examples",
    }
  }

  pub const fn env_var(self) -> &'static str {
    match self {
      Self::Profile => "TREETIME_PROFILE_DIR",
      Self::Runs => "TREETIME_RUNS_DIR",
      Self::Logs => "TREETIME_LOGS_DIR",
      Self::Examples => "TREETIME_EXAMPLES_DIR",
    }
  }
}

fn env_path(name: &str) -> Result<Option<PathBuf>, Report> {
  env_var_optional(name)?
    .filter(|value| !value.is_empty())
    .map(absolute_path)
    .transpose()
}
