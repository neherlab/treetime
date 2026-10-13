#[cfg(test)]
mod tests {
  use crate::app_paths::{AppFolderEnv, AppFolderPath, AppPaths};
  use crate::app_settings::settings::AppPathSettings;
  use helpers::unfixed;
  use pretty_assertions::assert_eq;
  use std::path::{Path, PathBuf};

  #[test]
  fn test_app_paths_default_to_folders_in_the_root() {
    let actual = AppPaths::resolve(
      Path::new("/home/alice/.config/treetime"),
      &AppFolderEnv::default(),
      &AppPathSettings::default(),
    );
    let expected = AppPaths {
      root: PathBuf::from("/home/alice/.config/treetime"),
      profile: unfixed("/home/alice/.config/treetime/profile"),
      runs: unfixed("/home/alice/.config/treetime/runs"),
      logs: unfixed("/home/alice/.config/treetime/logs"),
      examples: unfixed("/home/alice/.config/treetime/examples"),
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_app_paths_take_settings_relative_to_the_root_or_absolute() {
    let settings = AppPathSettings {
      runs: Some(PathBuf::from("/data/runs")),
      logs: Some(PathBuf::from("state/logs")),
      ..AppPathSettings::default()
    };
    let actual = AppPaths::resolve(Path::new("/app"), &AppFolderEnv::default(), &settings);
    assert_eq!(
      (
        unfixed("/data/runs"),
        unfixed("/app/state/logs"),
        unfixed("/app/profile")
      ),
      (actual.runs, actual.logs, actual.profile)
    );
  }

  #[test]
  fn test_app_paths_environment_wins_over_settings() {
    let env = AppFolderEnv {
      runs: Some(PathBuf::from("/scratch/runs")),
      examples: Some(PathBuf::from("/checkout/data")),
      ..AppFolderEnv::default()
    };
    let settings = AppPathSettings {
      runs: Some(PathBuf::from("/data/runs")),
      ..AppPathSettings::default()
    };
    let actual = AppPaths::resolve(Path::new("/app"), &env, &settings);
    let expected = (
      AppFolderPath {
        path: PathBuf::from("/scratch/runs"),
        fixed_by: Some("TREETIME_RUNS_DIR"),
      },
      AppFolderPath {
        path: PathBuf::from("/checkout/data"),
        fixed_by: Some("TREETIME_EXAMPLES_DIR"),
      },
    );
    assert_eq!(expected, (actual.runs, actual.examples));
  }

  #[test]
  fn test_app_paths_default_runs_ignores_settings_and_environment() {
    let env = AppFolderEnv {
      runs: Some(PathBuf::from("/scratch/runs")),
      ..AppFolderEnv::default()
    };
    let actual = AppPaths::resolve(Path::new("/app"), &env, &AppPathSettings::default());
    assert_eq!(PathBuf::from("/app/runs"), actual.default_runs());
  }

  mod helpers {
    use crate::app_paths::AppFolderPath;
    use std::path::PathBuf;

    pub(super) fn unfixed(path: &str) -> AppFolderPath {
      AppFolderPath {
        path: PathBuf::from(path),
        fixed_by: None,
      }
    }
  }
}
