#[cfg(test)]
mod tests {
  use crate::command::config_from_settings;
  use crate::job::JobId;
  use crate::runs::store::{NewRun, RunStore};
  use app_commands::command::AppCommand;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use tempfile::tempdir;

  #[test]
  fn test_store_create_recreates_a_removed_runs_directory() {
    let root = tempdir().unwrap();
    let runs_dir = root.path().join("app").join("runs");
    let store = RunStore::open(&runs_dir).unwrap();
    fs::remove_dir_all(root.path().join("app")).unwrap();

    let record = store
      .create(NewRun {
        id: JobId::random(),
        title: "prune".to_owned(),
        config: config_from_settings(AppCommand::Prune, &json!({ "tree": "t.nwk" })).unwrap(),
        inputs: vec![],
        source_run_id: None,
        example: None,
      })
      .unwrap();

    assert!(store.run_dir(&record.id).join("out").is_dir());
    assert_eq!(
      vec![record.id],
      store.list().unwrap().into_iter().map(|run| run.id).collect::<Vec<_>>()
    );
  }
}
