#[cfg(test)]
mod tests {
  use crate::app_settings::settings::{AnalysisSettings, AppPathSettings, AppSettings, UiSettings, UiTheme};
  use crate::app_settings::store::{AppSettingsStore, SETTINGS_JSON, SETTINGS_YAML};
  use helpers::{draft, read_text, runs_at};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::PathBuf;
  use tempfile::tempdir;
  use treetime_grid::MaxGridPoints;
  use treetime_utils::assert_error;

  #[test]
  fn test_store_reads_defaults_when_no_settings_file_exists() {
    let dir = tempdir().unwrap();
    let store = AppSettingsStore::open(dir.path()).unwrap();
    assert_eq!(
      (dir.path().join(SETTINGS_YAML), AppSettings::default()),
      (store.path().to_path_buf(), store.read().unwrap())
    );
  }

  #[test]
  fn test_store_writes_yaml_by_default() {
    let dir = tempdir().unwrap();
    let store = AppSettingsStore::open(dir.path()).unwrap();
    store
      .update(|settings| {
        settings.paths.runs = Some(PathBuf::from("/data/runs"));
        settings.ui.theme = Some(UiTheme::Dark);
        settings.ui.sidebar_width = Some(400);
      })
      .unwrap();
    let expected = indoc! {r#"
      paths:
        runs: "/data/runs"
      ui:
        theme: "dark"
        sidebar_width: 400
    "#};
    assert_eq!(expected, read_text(&dir.path().join(SETTINGS_YAML)));
  }

  #[test]
  fn test_store_writes_json_when_the_settings_file_is_json() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join(SETTINGS_JSON), "{}\n").unwrap();
    let store = AppSettingsStore::open(dir.path()).unwrap();
    store
      .update(|settings| settings.ui.theme = Some(UiTheme::Light))
      .unwrap();
    let expected = indoc! {r#"{
      "ui": {
        "theme": "light"
      }
    }
    "#};
    assert_eq!(
      (dir.path().join(SETTINGS_JSON), expected.to_owned()),
      (store.path().to_path_buf(), read_text(&dir.path().join(SETTINGS_JSON)))
    );
  }

  #[test]
  fn test_store_reads_back_a_draft_from_yaml() {
    let dir = tempdir().unwrap();
    let store = AppSettingsStore::open(dir.path()).unwrap();
    let expected = AppSettings {
      paths: AppPathSettings::default(),
      ui: UiSettings {
        theme: Some(UiTheme::System),
        sidebar_width: None,
        draft: Some(draft()),
      },
      analysis: AnalysisSettings {
        max_grid_points: Some(MaxGridPoints::new(250_000).unwrap()),
      },
    };
    store.update(|settings| *settings = expected.clone()).unwrap();
    assert_eq!(expected, AppSettingsStore::open(dir.path()).unwrap().read().unwrap());
  }

  #[test]
  fn test_store_reads_back_a_draft_from_json() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join(SETTINGS_JSON), "").unwrap();
    let store = AppSettingsStore::open(dir.path()).unwrap();
    let expected = AppSettings {
      paths: runs_at("/data/runs"),
      ui: UiSettings {
        theme: Some(UiTheme::Dark),
        sidebar_width: Some(300),
        draft: Some(draft()),
      },
      analysis: AnalysisSettings {
        max_grid_points: Some(MaxGridPoints::new(5_000).unwrap()),
      },
    };
    store.update(|settings| *settings = expected.clone()).unwrap();
    assert_eq!(expected, AppSettingsStore::open(dir.path()).unwrap().read().unwrap());
  }

  #[test]
  fn test_store_update_keeps_the_settings_it_does_not_change() {
    let dir = tempdir().unwrap();
    let store = AppSettingsStore::open(dir.path()).unwrap();
    store
      .update(|settings| settings.paths.runs = Some(PathBuf::from("/data/runs")))
      .unwrap();
    let updated = store
      .update(|settings| settings.ui.theme = Some(UiTheme::Dark))
      .unwrap();
    assert_eq!(
      AppSettings {
        paths: runs_at("/data/runs"),
        ui: UiSettings {
          theme: Some(UiTheme::Dark),
          ..UiSettings::default()
        },
        analysis: AnalysisSettings::default(),
      },
      updated
    );
  }

  #[test]
  fn test_store_reads_defaults_from_an_empty_file() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join(SETTINGS_YAML), "\n").unwrap();
    assert_eq!(
      AppSettings::default(),
      AppSettingsStore::open(dir.path()).unwrap().read().unwrap()
    );
  }

  #[test]
  fn test_store_refuses_both_a_yaml_and_a_json_file() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join(SETTINGS_YAML), "").unwrap();
    fs::write(dir.path().join(SETTINGS_JSON), "").unwrap();
    assert_error!(
      AppSettingsStore::open(dir.path()),
      format!(
        "both '{}' and '{}' exist; keep one of them",
        dir.path().join(SETTINGS_YAML).display(),
        dir.path().join(SETTINGS_JSON).display()
      )
    );
  }

  #[test]
  fn test_store_names_the_file_of_an_unknown_setting() {
    let dir = tempdir().unwrap();
    let path = dir.path().join(SETTINGS_YAML);
    fs::write(&path, "theme: dark\n").unwrap();
    let expected = indoc! {"
      error: line 1 column 1: unknown field `theme`, expected one of paths, ui
       --> <input>:1:1
        |
      1 | theme: dark
        | ^ unknown field `theme`, expected one of paths, ui"};
    assert_error!(
      AppSettingsStore::open(dir.path()).unwrap().read(),
      format!("When reading the settings file '{}': {expected}", path.display())
    );
  }

  #[test]
  fn test_store_names_the_file_of_a_duplicate_key() {
    let dir = tempdir().unwrap();
    let path = dir.path().join(SETTINGS_YAML);
    fs::write(&path, "ui:\n  theme: dark\nui:\n  theme: light\n").unwrap();
    let expected = indoc! {"
      error: line 3 column 1: duplicate mapping key: ui, set DuplicateKeyPolicy in Options if acceptable
       --> <input>:3:1
        |
      1 | ui:
      2 |   theme: dark
      3 | ui:
        | ^ duplicate mapping key: ui, set DuplicateKeyPolicy in Options if acceptable
      4 |   theme: light
        |"};
    assert_error!(
      AppSettingsStore::open(dir.path()).unwrap().read(),
      format!("When reading the settings file '{}': {expected}", path.display())
    );
  }

  #[test]
  fn test_store_names_the_file_of_malformed_json() {
    let dir = tempdir().unwrap();
    let path = dir.path().join(SETTINGS_JSON);
    fs::write(&path, "{\"ui\": {\"theme\": \"blue\"}}\n").unwrap();
    assert_error!(
      AppSettingsStore::open(dir.path()).unwrap().read(),
      format!(
        "When reading the settings file '{}': When parsing JSON: unknown variant `blue`, expected one of `system`, `light`, `dark` at line 1 column 23",
        path.display()
      )
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::sparse;
    use crate::app_settings::settings::{
      AppPathSettings, UiCodeFormat, UiDraft, UiDraftOrigin, UiDraftSource, UiSettingsView,
    };
    use crate::command::AppCommand;
    use crate::job::JobId;
    use maplit::btreemap;
    use serde_json::json;
    use std::fs;
    use std::path::{Path, PathBuf};
    use treetime_utils::o;

    pub(super) fn draft() -> UiDraft {
      UiDraft {
        command: AppCommand::Timetree,
        config: sparse(
          json!({ "tree": "data/flu/h3n2/20/tree.nwk", "clock_rate": 0.003, "dates": { "column": "date" } }),
        ),
        sources: btreemap! {
          o!("tree") => UiDraftSource { label: o!("tree.nwk"), origin: UiDraftOrigin::Dataset, size: Some(1234) },
          o!("alignment") => UiDraftSource { label: o!("aln: \"x\".fasta"), origin: UiDraftOrigin::Local, size: None },
        },
        from_run_id: Some(JobId::random()),
        upload_run_id: None,
        view: UiSettingsView::All,
        search: o!("clock: rate"),
        changed_only: true,
        code_format: UiCodeFormat::Yaml,
      }
    }

    pub(super) fn read_text(path: &Path) -> String {
      fs::read_to_string(path).unwrap()
    }

    pub(super) fn runs_at(path: &str) -> AppPathSettings {
      AppPathSettings {
        runs: Some(PathBuf::from(path)),
        ..AppPathSettings::default()
      }
    }
  }
}
