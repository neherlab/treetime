use crate::app_paths::AppPaths;
use crate::app_settings::settings::Workspace;
use crate::runs::errors::invalid;
use eyre::Report;
use std::fs;
use std::path::{Path, PathBuf};

pub fn active_workspace(paths: &AppPaths) -> Workspace {
  Workspace {
    path: paths.runs.path.clone(),
    default_path: paths.default_runs(),
    fixed_by: paths.runs.fixed_by.map(str::to_owned),
  }
}

pub fn prepare_workspace(paths: &AppPaths, path: &Path) -> Result<PathBuf, Report> {
  if let Some(variable) = paths.runs.fixed_by {
    return Err(invalid(format!(
      "the environment variable {variable} sets the runs folder; unset it to choose the folder in the app"
    )));
  }
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
