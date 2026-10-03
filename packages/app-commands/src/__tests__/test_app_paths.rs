#[cfg(test)]
mod tests {
  use crate::app_paths::{AppPaths, PlatformDirs};
  use pretty_assertions::assert_eq;
  use std::path::PathBuf;

  #[test]
  fn test_app_paths_in_dir_keeps_every_folder_below_the_root() {
    let actual = AppPaths::in_dir(&PathBuf::from("/checkout/tmp/app/treetime-dev"));
    let expected = AppPaths {
      profile_dir: PathBuf::from("/checkout/tmp/app/treetime-dev/profile"),
      settings_dir: PathBuf::from("/checkout/tmp/app/treetime-dev"),
      default_workspace: PathBuf::from("/checkout/tmp/app/treetime-dev/runs"),
      logs_dir: PathBuf::from("/checkout/tmp/app/treetime-dev/logs"),
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_app_paths_platform_names_every_folder_treetime() {
    let actual = AppPaths::from_platform_dirs(&PlatformDirs {
      config: PathBuf::from("/home/alice/.config"),
      data: PathBuf::from("/home/alice/.local/share"),
      logs: PathBuf::from("/home/alice/.local/state/treetime/logs"),
    });
    let expected = AppPaths {
      profile_dir: PathBuf::from("/home/alice/.config/treetime"),
      settings_dir: PathBuf::from("/home/alice/.config/treetime"),
      default_workspace: PathBuf::from("/home/alice/.local/share/treetime/runs"),
      logs_dir: PathBuf::from("/home/alice/.local/state/treetime/logs"),
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_app_paths_platform_shares_one_folder_when_config_and_data_coincide() {
    let actual = AppPaths::from_platform_dirs(&PlatformDirs {
      config: PathBuf::from("/Users/Alice/Library/Application Support"),
      data: PathBuf::from("/Users/Alice/Library/Application Support"),
      logs: PathBuf::from("/Users/Alice/Library/Logs/treetime"),
    });
    let expected = AppPaths {
      profile_dir: PathBuf::from("/Users/Alice/Library/Application Support/treetime"),
      settings_dir: PathBuf::from("/Users/Alice/Library/Application Support/treetime"),
      default_workspace: PathBuf::from("/Users/Alice/Library/Application Support/treetime/runs"),
      logs_dir: PathBuf::from("/Users/Alice/Library/Logs/treetime"),
    };
    assert_eq!(expected, actual);
  }

  #[cfg(all(unix, not(target_os = "macos")))]
  #[test]
  fn test_app_paths_platform_resolves_the_xdg_folders_of_the_user() {
    let actual = AppPaths::platform().unwrap();
    let config = dirs::config_local_dir().unwrap();
    let data = dirs::data_local_dir().unwrap();
    let state = dirs::state_dir().unwrap();
    let expected = AppPaths {
      profile_dir: config.join("treetime"),
      settings_dir: config.join("treetime"),
      default_workspace: data.join("treetime").join("runs"),
      logs_dir: state.join("treetime").join("logs"),
    };
    assert_eq!(expected, actual);
  }
}
