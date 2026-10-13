use crate::app_paths::{AppFolderPath, AppPaths};
use crate::app_settings::settings::Workspace;
use crate::runs::errors::invalid;
use crate::runs::manager::RunManager;
use eyre::Report;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use treetime_utils::error::report_to_string;
use treetime_utils::make_report;

pub fn active_workspace(paths: &AppPaths, error: Option<&str>) -> Workspace {
  Workspace {
    path: paths.runs.path.clone(),
    default_path: paths.default_runs(),
    fixed_by: paths.runs.fixed_by.map(str::to_owned),
    error: error.map(str::to_owned),
  }
}

pub struct OpenedRuns {
  pub runs: Arc<RunManager>,
  pub error: Option<String>,
}

pub fn open_runs_folder(paths: &mut AppPaths, named_in: Option<&Path>) -> Result<OpenedRuns, Report> {
  let source = |folder: &AppFolderPath| match (folder.fixed_by, named_in) {
    (Some(variable), _) => format!(", set by {variable}"),
    (None, Some(settings_file)) => format!(", named in '{}'", settings_file.display()),
    (None, None) => String::new(),
  };
  let report = match RunManager::open(&paths.runs.path) {
    Ok(runs) => return Ok(OpenedRuns { runs, error: None }),
    Err(report) => report,
  };
  let failed = format!("the runs folder '{}'{}", paths.runs.path.display(), source(&paths.runs));
  if paths.runs.fixed_by.is_some() || named_in.is_none() {
    return Err(report.wrap_err(format!("When opening {failed}")));
  }
  let error = format!("{failed} cannot be opened: {}", report_to_string(&report));
  let default = paths.default_runs();
  let runs = RunManager::open(&default).map_err(|fallback| {
    make_report!(
      "{error}; the default runs folder '{}' cannot be opened either: {}",
      default.display(),
      report_to_string(&fallback)
    )
  })?;
  paths.runs = AppFolderPath {
    path: default,
    fixed_by: None,
  };
  Ok(OpenedRuns {
    runs,
    error: Some(format!("{error}. The default runs folder is in use.")),
  })
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
