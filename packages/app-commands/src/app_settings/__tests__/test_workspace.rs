#[cfg(test)]
mod tests {
  use crate::app_paths::{AppFolderEnv, AppPaths};
  use crate::app_settings::settings::{AppPathSettings, Workspace};
  use crate::app_settings::workspace::{active_workspace, prepare_workspace};
  use helpers::default_paths;
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_workspace_active_is_the_runs_folder_of_the_settings() {
    let settings = AppPathSettings {
      runs: Some(PathBuf::from("/data/runs")),
      ..AppPathSettings::default()
    };
    let paths = AppPaths::resolve(Path::new("/app"), &AppFolderEnv::default(), &settings);
    let expected = Workspace {
      path: PathBuf::from("/data/runs"),
      default_path: PathBuf::from("/app/runs"),
      fixed_by: None,
    };
    assert_eq!(expected, active_workspace(&paths));
  }

  #[test]
  fn test_workspace_active_names_the_variable_that_fixes_it() {
    let env = AppFolderEnv {
      runs: Some(PathBuf::from("/scratch/runs")),
      ..AppFolderEnv::default()
    };
    let paths = AppPaths::resolve(Path::new("/app"), &env, &AppPathSettings::default());
    let expected = Workspace {
      path: PathBuf::from("/scratch/runs"),
      default_path: PathBuf::from("/app/runs"),
      fixed_by: Some("TREETIME_RUNS_DIR".to_owned()),
    };
    assert_eq!(expected, active_workspace(&paths));
  }

  #[test]
  fn test_workspace_prepare_creates_a_missing_folder() {
    let dir = tempdir().unwrap();
    let folder = dir.path().join("a").join("runs");
    let prepared = prepare_workspace(&default_paths(), &folder).unwrap();
    assert_eq!((folder.clone(), true), (prepared, folder.is_dir()));
  }

  #[test]
  fn test_workspace_prepare_refuses_a_relative_path() {
    assert_error!(
      prepare_workspace(&default_paths(), Path::new("runs")),
      "the runs folder must be an absolute path, got 'runs'"
    );
  }

  #[test]
  fn test_workspace_prepare_refuses_a_folder_fixed_by_the_environment() {
    let env = AppFolderEnv {
      runs: Some(PathBuf::from("/scratch/runs")),
      ..AppFolderEnv::default()
    };
    let paths = AppPaths::resolve(Path::new("/app"), &env, &AppPathSettings::default());
    assert_error!(
      prepare_workspace(&paths, Path::new("/data/runs")),
      "the environment variable TREETIME_RUNS_DIR sets the runs folder; unset it to choose the folder in the app"
    );
  }

  #[cfg(target_os = "linux")]
  #[test]
  fn test_workspace_prepare_refuses_a_file() {
    let dir = tempdir().unwrap();
    let file = dir.path().join("runs");
    fs::write(&file, "").unwrap();
    assert_error!(
      prepare_workspace(&default_paths(), &file),
      format!(
        "the runs folder '{}' cannot be created: File exists (os error 17)",
        file.display()
      )
    );
  }

  mod helpers {
    use crate::app_paths::{AppFolderEnv, AppPaths};
    use crate::app_settings::settings::AppPathSettings;
    use std::path::Path;

    pub(super) fn default_paths() -> AppPaths {
      AppPaths::resolve(Path::new("/app"), &AppFolderEnv::default(), &AppPathSettings::default())
    }
  }
}
