#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{sparse, unconfined};
  use crate::config::code::config_code;
  use crate::results::run_results::run_results;
  use crate::runs::app_events::AppChange;
  use crate::runs::manager::RunManager;
  use crate::runs::record::{CreateRunRequest, RunStatus, StartRunRequest};
  use app_commands::command::AppCommand;
  use helpers::{ZIKA, run_with_uploads, settings, start_from};
  use itertools::Itertools;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use std::sync::Arc;
  use tempfile::tempdir;
  use treetime_utils::{assert_error, o};

  #[test]
  fn test_run_storage_record_holds_paths_relative_to_the_run_folder() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));

    let settings = settings(&record);
    assert_eq!(
      (
        RunStatus::Ok,
        json!("inputs/tree.nwk"),
        json!("inputs/metadata.tsv"),
        json!("out"),
        vec![o!("inputs/metadata.tsv"), o!("inputs/tree.nwk")],
      ),
      (
        record.status,
        settings["tree"].clone(),
        settings["metadata"].clone(),
        settings["output_all"].clone(),
        record
          .inputs
          .iter()
          .map(|input| input.path.to_string_lossy().into_owned())
          .sorted()
          .collect_vec(),
      )
    );
  }

  #[test]
  fn test_run_storage_relative_paths_survive_a_move_of_the_runs_root() {
    let root = tempdir().unwrap();
    let before = root.path().join("before");
    let id = {
      let runs = RunManager::open(&before).unwrap();
      run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({})).id
    };
    let after = root.path().join("after");
    fs::rename(&before, &after).unwrap();
    let runs = RunManager::open(&after).unwrap();

    let rerun = start_from(&runs, &runs.settings(&id).unwrap(), Some(&id));
    let results = run_results(&runs, &id).unwrap();

    assert_eq!(
      (true, RunStatus::Ok, true),
      (
        runs.settings(&id).unwrap().0["tree"]
          .as_str()
          .unwrap()
          .starts_with(&after.canonicalize().unwrap().to_string_lossy().into_owned()),
        rerun.status,
        results.problems.is_empty(),
      )
    );
  }

  #[test]
  fn test_run_storage_variant_keeps_working_after_the_delete_of_its_source() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let source = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let variant = start_from(&runs, &runs.settings(&source.id).unwrap(), Some(&source.id));

    runs.delete(&source.id).unwrap();
    let rerun = start_from(&runs, &runs.settings(&variant.id).unwrap(), Some(&variant.id));

    assert_eq!(
      (
        (RunStatus::Ok, Some(source.id.clone()), source.group_id.clone()),
        (RunStatus::Ok, source.group_id),
        true,
      ),
      (
        (variant.status, variant.source_run_id.clone(), variant.group_id.clone()),
        (rerun.status, rerun.group_id),
        run_results(&runs, &variant.id).unwrap().problems.is_empty(),
      )
    );
  }

  #[test]
  fn test_run_storage_delete_refuses_a_computing_run() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = runs
      .create(&CreateRunRequest {
        command: AppCommand::Clock,
        config: sparse(json!({})),
        defer_start: true,
        source_run_id: None,
      })
      .unwrap();
    let started = runs
      .start(&created.id, &StartRunRequest::default(), unconfined())
      .unwrap();

    assert_error!(
      runs.delete(&created.id),
      format!(
        "run `{}` is computing; cancel it before deleting it",
        created.id.as_str()
      )
    );
    assert_eq!(RunStatus::Running, runs.get(&created.id).unwrap().status);
    drop(started);
  }

  #[test]
  fn test_run_storage_delete_emits_one_event_and_a_second_delete_finds_no_run() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let changes = Arc::new(Mutex::new(vec![]));
    let seen = Arc::clone(&changes);
    runs.app_events().subscribe(
      None,
      Box::new(move |event| {
        seen.lock().push(event.change.clone());
        true
      }),
    );

    runs.delete(&record.id).unwrap();

    assert_error!(
      runs.delete(&record.id),
      format!("no run with id `{}`", record.id.as_str())
    );
    let changes = changes.lock();
    let [AppChange::RunRemoved { id }] = changes.as_slice() else {
      panic!("expected one run-removed change, got {changes:?}");
    };
    assert_eq!(
      (&record.id, false, 0),
      (
        id,
        runs.store().run_dir(&record.id).exists(),
        runs.input_store().keys().unwrap().len()
      )
    );
  }

  #[test]
  fn test_run_storage_delete_keeps_the_inputs_that_another_run_references() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let source = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let variant = start_from(&runs, &runs.settings(&source.id).unwrap(), Some(&source.id));
    let shared_before = runs.get(&variant.id).unwrap().size;

    runs.delete(&source.id).unwrap();
    let alone = runs.get(&variant.id).unwrap().size;
    let tree = runs.settings(&variant.id).unwrap().0["tree"].clone();
    runs.delete(&variant.id).unwrap();

    assert_eq!(
      (true, false, 0),
      (
        alone > shared_before,
        fs::metadata(tree.as_str().unwrap()).is_ok(),
        runs.input_store().keys().unwrap().len(),
      )
    );
  }

  #[test]
  fn test_run_storage_run_without_a_source_joins_the_group_of_an_earlier_run_with_the_same_inputs() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let first = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    let same_inputs = start_from(&runs, &runs.settings(&first.id).unwrap(), None);
    let other_inputs = run_with_uploads(&runs, AppCommand::Clock, "ebola/20", &json!({}));

    assert_eq!(
      (first.id.clone(), first.id.clone(), other_inputs.id.clone()),
      (first.group_id, same_inputs.group_id, other_inputs.group_id)
    );
  }

  #[test]
  fn test_run_storage_records_the_drawn_seed_and_the_config_code_carries_it() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(
      &runs,
      AppCommand::Timetree,
      ZIKA,
      &json!({ "max_iter": 1, "sequence_length": 10807 }),
    );
    assert_eq!(
      (RunStatus::Ok, None),
      (record.status, record.error.as_ref().map(|error| &error.message))
    );
    let seed = record
      .seed
      .expect("a time tree with polytomy resolution records its seed");
    let code = config_code(AppCommand::Timetree, &record.config.settings().unwrap()).unwrap();

    assert_eq!(
      (json!(seed), false, true, true),
      (
        settings(&record)["seed"].clone(),
        record.changed_settings.contains(&"seed".to_owned()),
        code.command_line_text.contains(&format!("--seed {seed}")),
        code.command_line_text.contains("--tree inputs/tree.nwk"),
      )
    );
  }

  #[test]
  fn test_run_storage_records_no_seed_for_a_run_without_a_random_step() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let record = run_with_uploads(&runs, AppCommand::Clock, ZIKA, &json!({}));
    assert_eq!(
      (RunStatus::Ok, None, None),
      (record.status, record.seed, settings(&record).get("seed"))
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::{sparse, unconfined};
    use crate::job::JobId;
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunRecord, StartRunRequest};
    use app_commands::command::AppCommand;
    use app_commands::json_value::SparseConfig;
    use serde_json::{Map, Value};
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::Arc;

    pub(super) const ZIKA: &str = "zika/20";

    pub(super) fn run_with_uploads(
      runs: &Arc<RunManager>,
      command: AppCommand,
      dataset: &str,
      extra: &Value,
    ) -> RunRecord {
      let created = runs
        .create(&CreateRunRequest {
          command,
          config: sparse(Value::Object(Map::new())),
          defer_start: true,
          source_run_id: None,
        })
        .unwrap();
      let mut config = extra.as_object().cloned().unwrap_or_default();
      for (key, name) in [("tree", "tree.nwk"), ("metadata", "metadata.tsv")] {
        let path = upload(runs, &created.id, dataset, name);
        config.insert(key.to_owned(), Value::from(path.to_string_lossy().into_owned()));
      }
      let request = StartRunRequest {
        command: None,
        config: Some(SparseConfig(config)),
        source_run_id: None,
      };
      runs.start(&created.id, &request, unconfined()).unwrap().run();
      runs.get(&created.id).unwrap()
    }

    fn upload(runs: &RunManager, id: &JobId, dataset: &str, name: &str) -> PathBuf {
      let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data").join(dataset);
      let bytes = fs::read(data.join(name)).unwrap();
      runs.upload_input(id, name, &[], &mut &bytes[..], 1 << 24).unwrap().path
    }

    pub(super) fn start_from(runs: &Arc<RunManager>, settings: &SparseConfig, source: Option<&JobId>) -> RunRecord {
      let created = runs
        .create(&CreateRunRequest {
          command: AppCommand::Clock,
          config: settings.clone(),
          defer_start: true,
          source_run_id: source.cloned(),
        })
        .unwrap();
      runs
        .start(&created.id, &StartRunRequest::default(), unconfined())
        .unwrap()
        .run();
      runs.get(&created.id).unwrap()
    }

    pub(super) fn settings(record: &RunRecord) -> Map<String, Value> {
      record.config.settings().unwrap()
    }
  }
}
