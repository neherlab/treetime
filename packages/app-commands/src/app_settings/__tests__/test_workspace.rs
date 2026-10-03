#[cfg(test)]
mod tests {
  use crate::app_paths::AppPaths;
  use crate::app_settings::settings::{AppSettings, Workspace};
  use crate::app_settings::workspace::{active_workspace, prepare_workspace};
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::{Path, PathBuf};
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_workspace_active_is_the_default_folder_without_a_setting() {
    let paths = AppPaths::in_dir(Path::new("/app"));
    let expected = Workspace {
      path: PathBuf::from("/app/runs"),
      default_path: PathBuf::from("/app/runs"),
    };
    assert_eq!(expected, active_workspace(&AppSettings::default(), &paths));
  }

  #[test]
  fn test_workspace_active_is_the_folder_the_settings_name() {
    let paths = AppPaths::in_dir(Path::new("/app"));
    let settings = AppSettings {
      workspace: Some(PathBuf::from("/data/runs")),
      ..AppSettings::default()
    };
    let expected = Workspace {
      path: PathBuf::from("/data/runs"),
      default_path: PathBuf::from("/app/runs"),
    };
    assert_eq!(expected, active_workspace(&settings, &paths));
  }

  #[test]
  fn test_workspace_prepare_creates_a_missing_folder() {
    let dir = tempdir().unwrap();
    let folder = dir.path().join("a").join("runs");
    let prepared = prepare_workspace(&folder).unwrap();
    assert_eq!((folder.clone(), true), (prepared, folder.is_dir()));
  }

  #[test]
  fn test_workspace_prepare_refuses_a_relative_path() {
    assert_error!(
      prepare_workspace(Path::new("runs")),
      "the runs folder must be an absolute path, got 'runs'"
    );
  }

  #[cfg(target_os = "linux")]
  #[test]
  fn test_workspace_prepare_refuses_a_file() {
    let dir = tempdir().unwrap();
    let file = dir.path().join("runs");
    fs::write(&file, "").unwrap();
    assert_error!(
      prepare_workspace(&file),
      format!(
        "the runs folder '{}' cannot be created: File exists (os error 17)",
        file.display()
      )
    );
  }
}
