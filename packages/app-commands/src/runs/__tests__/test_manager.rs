#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::job::{JobEvent, TerminalEvent};
  use crate::runs::events::{RunEvent, read_events};
  use crate::runs::manager::RunManager;
  use crate::runs::record::{RunRecord, RunStatus, UpdateRunRequest};
  use crate::runs::store::RunStore;
  use helpers::{accept, clock_config, create, event_types, terminal_status, timetree_config};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use std::sync::Arc;
  use std::thread;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_manager_run_records_its_lifecycle_and_outputs() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    assert_eq!(RunStatus::Created, created.status);

    let terminal = runs.start(&created.id, None, accept()).unwrap().run();
    assert_eq!("ok", terminal_status(&terminal));

    let record = runs.get(&created.id).unwrap();
    let out_dir = runs.store().out_dir(&created.id);
    assert_eq!(
      (RunStatus::Ok, json!(out_dir), true, true, vec!["started", "terminal"],),
      (
        record.status,
        record.config["output_all"].clone(),
        record
          .output_files
          .iter()
          .all(|file| out_dir.join(&file.path).is_file()),
        record.headline.contains_key("clock_rate") && record.headline.contains_key("r_squared"),
        event_types(&read_events(&runs.store().events_path(&created.id), 0).unwrap())
          .into_iter()
          .filter(|kind| *kind == "started" || *kind == "terminal")
          .collect::<Vec<_>>(),
      )
    );
    assert!(record.started_at.is_some() && record.finished_at.is_some() && record.duration_seconds.is_some());
    assert_eq!(2, record.inputs.len(), "the tree and the metadata are recorded");
  }

  #[test]
  fn test_manager_run_adds_the_auspice_output_and_keeps_the_chosen_selection() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let mut config = clock_config();
    config["output_selection"] = json!(["ClockModel"]);
    config["output_all"] = json!("/elsewhere");
    config["output_clock_model"] = json!("/elsewhere/model.json");
    let created = create(&runs, AppCommand::Clock, config);
    runs.start(&created.id, None, accept()).unwrap().run();
    let record = runs.get(&created.id).unwrap();
    let files: Vec<&str> = record
      .output_files
      .iter()
      .map(|file| file.path.to_str().unwrap())
      .collect();
    assert_eq!(
      (
        json!(["ClockModel", "Auspice"]),
        json!(null),
        vec!["clock.auspice.json", "clock.clock-model.json"]
      ),
      (
        record.config["output_selection"].clone(),
        record.config["output_clock_model"].clone(),
        files
      )
    );
    assert_eq!(vec!["output_selection".to_owned()], record.changed_settings);
  }

  #[test]
  fn test_manager_lists_newest_first_and_counts_active_runs() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let first = create(&runs, AppCommand::Clock, clock_config());
    let second = create(&runs, AppCommand::Clock, clock_config());
    let started = runs.start(&second.id, None, accept()).unwrap();
    let list = runs.list().unwrap();
    assert_eq!(
      (vec![second.id, first.id], 1),
      (
        list.runs.iter().map(|run| run.id.clone()).collect::<Vec<_>>(),
        list.active_runs
      )
    );
    started.run();
    assert_eq!(0, runs.list().unwrap().active_runs);
  }

  #[test]
  fn test_manager_cancel_before_start_ends_the_run_with_a_cancelled_event() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    assert!(runs.cancel(&created.id).unwrap());
    let events = read_events(&runs.store().events_path(&created.id), 0).unwrap();
    assert_eq!(
      (RunStatus::Cancelled, vec!["terminal"], false),
      (
        runs.get(&created.id).unwrap().status,
        event_types(&events),
        runs.cancel(&created.id).unwrap()
      )
    );
    assert_error!(
      runs.start(&created.id, None, accept()),
      format!("run `{}` has already started", created.id.as_str())
    );
  }

  #[test]
  fn test_manager_two_concurrent_runs_cancel_independently() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let first = create(&runs, AppCommand::Timetree, timetree_config());
    let second = create(&runs, AppCommand::Timetree, timetree_config());
    let first_started = runs.start(&first.id, None, accept()).unwrap();
    let second_started = runs.start(&second.id, None, accept()).unwrap();

    let cancelled = Arc::new(Mutex::new(false));
    let manager = Arc::clone(&runs);
    let first_id = first.id.clone();
    let flag = Arc::clone(&cancelled);
    runs
      .subscribe(
        &first.id,
        0,
        Box::new(move |event: &RunEvent| {
          if matches!(event.event, JobEvent::Progress(_)) && !*flag.lock() {
            *flag.lock() = true;
            manager.cancel(&first_id).unwrap();
          }
          true
        }),
      )
      .unwrap();

    let (first_terminal, second_terminal) = thread::scope(|scope| {
      let first = scope.spawn(|| first_started.run());
      let second = scope.spawn(|| second_started.run());
      (first.join().unwrap(), second.join().unwrap())
    });
    let progress_of = |id| {
      read_events(&runs.store().events_path(id), 0)
        .unwrap()
        .iter()
        .filter(|event| matches!(event.event, JobEvent::Progress(_)))
        .count()
    };
    assert_eq!(
      ("cancelled", "ok", RunStatus::Cancelled, RunStatus::Ok),
      (
        terminal_status(&first_terminal),
        terminal_status(&second_terminal),
        runs.get(&first.id).unwrap().status,
        runs.get(&second.id).unwrap().status,
      )
    );
    assert!(
      progress_of(&second.id) > progress_of(&first.id),
      "each run has its own progress"
    );
  }

  #[test]
  fn test_manager_marks_runs_of_a_stopped_process_as_interrupted() {
    let root = tempdir().unwrap();
    let id = {
      let runs = RunManager::open(root.path()).unwrap();
      let created = create(&runs, AppCommand::Clock, clock_config());
      let started = runs.start(&created.id, None, accept()).unwrap();
      drop(started);
      created.id
    };
    let runs = RunManager::open(root.path()).unwrap();
    let events = read_events(&runs.store().events_path(&id), 0).unwrap();
    assert_eq!(
      (RunStatus::Interrupted, vec!["terminal"]),
      (runs.get(&id).unwrap().status, event_types(&events))
    );
    assert!(matches!(
      events.last().unwrap().event,
      JobEvent::Terminal(TerminalEvent::Interrupted { .. })
    ));
  }

  #[test]
  fn test_manager_delete_restore_and_purge() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    runs.delete(&created.id).unwrap();
    assert_error!(
      runs.get(&created.id),
      format!("no run with id `{}`", created.id.as_str())
    );
    assert!(runs.list().unwrap().runs.is_empty());
    assert_eq!(created.id, runs.restore(&created.id).unwrap().id);
    runs.delete(&created.id).unwrap();
    runs.purge(&created.id).unwrap();
    assert_error!(
      runs.restore(&created.id),
      format!("no deleted run with id `{}`", created.id.as_str())
    );
  }

  #[test]
  fn test_manager_refuses_to_delete_a_running_run() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    let started = runs.start(&created.id, None, accept()).unwrap();
    assert_error!(
      runs.delete(&created.id),
      format!("run `{}` is running; cancel it before deleting it", created.id.as_str())
    );
    started.run();
    runs.delete(&created.id).unwrap();
  }

  #[test]
  fn test_manager_rename_and_pin() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    let summary = runs
      .update(
        &created.id,
        UpdateRunRequest {
          title: Some("dengue".to_owned()),
          pinned: Some(true),
        },
      )
      .unwrap();
    assert_eq!(("dengue", true), (summary.title.as_str(), summary.pinned));
    assert_error!(
      runs.update(
        &created.id,
        UpdateRunRequest {
          title: Some(" ".to_owned()),
          pinned: None,
        },
      ),
      "a run title must not be empty"
    );
  }

  #[test]
  fn test_manager_upload_limit_counts_every_input_of_the_run() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let uploaded = runs
      .upload_input(&created.id, "tree.nwk", &mut &[b'x'; 700][..], 1024)
      .unwrap();
    assert_eq!(
      (700, runs.store().inputs_dir(&created.id).join("tree.nwk")),
      (uploaded.size, uploaded.path)
    );
    assert_error!(
      runs.upload_input(&created.id, "metadata.tsv", &mut &[b'x'; 400][..], 1024),
      "the inputs of a run are limited to 1 KiB in total; `metadata.tsv` does not fit"
    );
    assert!(!runs.store().inputs_dir(&created.id).join("metadata.tsv").exists());
    assert_eq!(
      1000,
      runs
        .upload_input(&created.id, "tree.nwk", &mut &[b'x'; 1000][..], 1024)
        .unwrap()
        .size,
      "replacing a file counts only its new size"
    );
  }

  #[test]
  fn test_manager_upload_rejects_names_with_directories() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    for name in ["../x.nwk", "a/b.nwk", ".hidden", ".."] {
      assert_error!(
        runs.upload_input(&created.id, name, &mut &b"x"[..], 1024),
        format!("upload name `{name}` must be a plain file name without directories")
      );
    }
  }

  #[test]
  fn test_manager_resolves_output_files_inside_the_output_folder_only() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    runs.start(&created.id, None, accept()).unwrap().run();
    fs::write(root.path().join("outside.txt"), "secret").unwrap();
    runs.file_path(&created.id, "clock.clock-model.json").unwrap();
    for escape in ["../run.json", "/etc/passwd", "../../outside.txt", ""] {
      assert_error!(
        runs.file_path(&created.id, escape),
        format!("file path `{escape}` must name a file inside the run's output folder")
      );
    }
    assert!(runs.zip(&created.id).unwrap().starts_with(b"PK"));
  }

  #[test]
  fn test_manager_store_writes_run_json_atomically() {
    let root = tempdir().unwrap();
    let store = RunStore::open(root.path()).unwrap();
    let record = store
      .create(AppCommand::Prune, json!({ "tree": "t.nwk" }), None)
      .unwrap();
    let dir = store.run_dir(&record.id);
    let names: Vec<String> = fs::read_dir(&dir)
      .unwrap()
      .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
      .collect::<std::collections::BTreeSet<_>>()
      .into_iter()
      .collect();
    let reread: RunRecord = serde_json::from_str(&fs::read_to_string(dir.join("run.json")).unwrap()).unwrap();
    assert_eq!(
      (vec!["inputs", "out", "run.json"], "prune"),
      (
        names.iter().map(String::as_str).collect::<Vec<_>>(),
        reread.title.as_str()
      )
    );
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::job::TerminalEvent;
    use crate::runs::events::RunEvent;
    use crate::runs::manager::{ConfigHook, RunManager};
    use crate::runs::record::{CreateRunRequest, RunRecord};
    use serde_json::{Value, json};
    use std::path::Path;

    pub(super) fn accept() -> ConfigHook {
      Box::new(|_config: &mut Value| Ok(()))
    }

    pub(super) fn create(runs: &RunManager, command: AppCommand, config: Value) -> RunRecord {
      runs
        .create(CreateRunRequest {
          command,
          config,
          title: None,
          defer_start: true,
        })
        .unwrap()
    }

    pub(super) fn clock_config() -> Value {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({ "tree": zika.join("tree.nwk"), "metadata": zika.join("metadata.tsv") })
    }

    pub(super) fn timetree_config() -> Value {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({
        "tree": zika.join("tree.nwk"),
        "metadata": zika.join("metadata.tsv"),
        "alignment": [zika.join("aln.fasta.xz")],
        "max_iter": 2,
        "seed": 7,
      })
    }

    pub(super) fn terminal_status(terminal: &TerminalEvent) -> &'static str {
      match terminal {
        TerminalEvent::Ok { .. } => "ok",
        TerminalEvent::Error { .. } => "error",
        TerminalEvent::Cancelled { .. } => "cancelled",
        TerminalEvent::Interrupted { .. } => "interrupted",
      }
    }

    pub(super) fn event_types(events: &[RunEvent]) -> Vec<&'static str> {
      events
        .iter()
        .map(
          |event| match serde_json::to_value(event).unwrap()["type"].as_str().unwrap() {
            "started" => "started",
            "progress" => "progress",
            "log" => "log",
            "iteration" => "iteration",
            "terminal" => "terminal",
            other => panic!("unknown event type {other}"),
          },
        )
        .collect()
    }
  }
}
