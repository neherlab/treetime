#[cfg(test)]
mod tests {
  use crate::command::config_from_settings;
  use crate::runs::record::RunRecord;
  use crate::runs::store::{RunStore, default_title};
  use app_commands::command::AppCommand;
  use chrono::{FixedOffset, Local, TimeZone};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::json;
  use std::collections::BTreeSet;
  use std::fs;
  use tempfile::tempdir;
  use treetime_utils::io::json::json_read_str;

  #[rustfmt::skip]
  #[rstest]
  #[case::utc(         0,         (2026, 9, 28, 14, 3, 59), "Run 2026-09-28 14:03")]
  #[case::east_of_utc( 2 * 3600,  (2026, 1,  2,  0, 0,  0), "Run 2026-01-02 00:00")]
  #[case::west_of_utc(-5 * 3600,  (2025, 12, 31, 23, 59, 0), "Run 2025-12-31 23:59")]
  #[trace]
  fn test_store_default_title_is_the_local_creation_minute(
    #[case] offset_seconds: i32,
    #[case] (year, month, day, hour, minute, second): (i32, u32, u32, u32, u32, u32),
    #[case] expected: &str,
  ) {
    let zone = FixedOffset::east_opt(offset_seconds).unwrap();
    let created_at = zone.with_ymd_and_hms(year, month, day, hour, minute, second).unwrap();
    assert_eq!(expected, default_title(&created_at));
  }

  #[test]
  fn test_store_create_recreates_a_removed_runs_directory() {
    let root = tempdir().unwrap();
    let runs_dir = root.path().join("app").join("runs");
    let store = RunStore::open(&runs_dir).unwrap();
    fs::remove_dir_all(root.path().join("app")).unwrap();

    let record = store
      .create(config_from_settings(AppCommand::Prune, &json!({ "tree": "t.nwk" })).unwrap())
      .unwrap();

    assert!(store.run_dir(&record.id).join("inputs").is_dir());
    assert_eq!(
      vec![record.id],
      store.list().unwrap().into_iter().map(|run| run.id).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_store_create_writes_run_json_without_leftover_files() {
    let root = tempdir().unwrap();
    let store = RunStore::open(root.path()).unwrap();
    let record = store
      .create(config_from_settings(AppCommand::Prune, &json!({ "tree": "t.nwk" })).unwrap())
      .unwrap();
    let dir = store.run_dir(&record.id);
    let names: Vec<String> = fs::read_dir(&dir)
      .unwrap()
      .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
      .collect::<BTreeSet<_>>()
      .into_iter()
      .collect();
    let reread: RunRecord = json_read_str(fs::read_to_string(dir.join("run.json")).unwrap()).unwrap();
    assert_eq!(
      (
        vec!["inputs", "out", "run.json"],
        default_title(&record.created_at.with_timezone(&Local))
      ),
      (names.iter().map(String::as_str).collect::<Vec<_>>(), reread.title)
    );
  }
}
