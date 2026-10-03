#[cfg(target_os = "macos")]
use dirs::home_dir;
#[cfg(not(any(target_os = "macos", windows)))]
use dirs::state_dir;
use dirs::{config_local_dir, data_local_dir};
use eyre::Report;
use std::path::{Path, PathBuf};
use treetime_utils::env::env_var_optional;
use treetime_utils::io::fs::absolute_path;
use treetime_utils::make_report;

pub const APP_DIR_ENV: &str = "TREETIME_APP_DIR";

const APP_NAME: &str = "treetime";

const PROFILE_DIR: &str = "profile";

const RUNS_DIR: &str = "runs";

const LOGS_DIR: &str = "logs";

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AppPaths {
  pub profile_dir: PathBuf,
  pub settings_dir: PathBuf,
  pub default_workspace: PathBuf,
  pub logs_dir: PathBuf,
}

impl AppPaths {
  pub fn from_env() -> Result<Self, Report> {
    match env_var_optional(APP_DIR_ENV)? {
      Some(dir) if !dir.is_empty() => Ok(Self::in_dir(&absolute_path(dir)?)),
      _ => Self::platform(),
    }
  }

  pub fn in_dir(root: &Path) -> Self {
    Self {
      profile_dir: root.join(PROFILE_DIR),
      settings_dir: root.to_path_buf(),
      default_workspace: root.join(RUNS_DIR),
      logs_dir: root.join(LOGS_DIR),
    }
  }

  pub fn platform() -> Result<Self, Report> {
    Ok(Self::from_platform_dirs(&PlatformDirs {
      config: user_dir(config_local_dir(), "configuration")?,
      data: user_dir(data_local_dir(), "data")?,
      logs: platform_logs_dir()?,
    }))
  }

  pub fn from_platform_dirs(platform: &PlatformDirs) -> Self {
    let config = platform.config.join(APP_NAME);
    Self {
      profile_dir: config.clone(),
      settings_dir: config,
      default_workspace: platform.data.join(APP_NAME).join(RUNS_DIR),
      logs_dir: platform.logs.clone(),
    }
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlatformDirs {
  pub config: PathBuf,
  pub data: PathBuf,
  pub logs: PathBuf,
}

fn user_dir(dir: Option<PathBuf>, kind: &str) -> Result<PathBuf, Report> {
  dir.ok_or_else(|| make_report!("the {kind} directory of the user could not be determined"))
}

#[cfg(target_os = "macos")]
fn platform_logs_dir() -> Result<PathBuf, Report> {
  Ok(
    user_dir(home_dir(), "home")?
      .join("Library")
      .join("Logs")
      .join(APP_NAME),
  )
}

#[cfg(windows)]
fn platform_logs_dir() -> Result<PathBuf, Report> {
  Ok(user_dir(data_local_dir(), "data")?.join(APP_NAME).join(LOGS_DIR))
}

#[cfg(not(any(target_os = "macos", windows)))]
fn platform_logs_dir() -> Result<PathBuf, Report> {
  Ok(user_dir(state_dir(), "state")?.join(APP_NAME).join(LOGS_DIR))
}
