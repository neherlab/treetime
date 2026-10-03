use crate::app_paths::AppPaths;
use crate::app_settings::settings::{AppSettings, Workspace};
use crate::runs::errors::invalid;
use eyre::Report;
use std::fs;
use std::path::{Path, PathBuf};

pub fn active_workspace(settings: &AppSettings, paths: &AppPaths) -> Workspace {
  Workspace {
    path: settings
      .workspace
      .clone()
      .unwrap_or_else(|| paths.default_workspace.clone()),
    default_path: paths.default_workspace.clone(),
  }
}

pub fn prepare_workspace(path: &Path) -> Result<PathBuf, Report> {
  if !path.is_absolute() {
    return Err(invalid(format!(
      "the runs folder must be an absolute path, got '{}'",
      path.display()
    )));
  }
  fs::create_dir_all(path)
    .map_err(|err| invalid(format!("the runs folder '{}' cannot be created: {err}", path.display())))?;
  Ok(path.to_path_buf())
}
