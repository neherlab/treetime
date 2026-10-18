#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::sparse;
  use crate::command::config_from_settings;
  use crate::job::{JobEvent, JobId, TerminalEvent};
  use crate::runs::events::{RunEvent, read_events};
  use crate::runs::manager::RunManager;
  use crate::runs::record::{RunRecord, RunStatus, StartRunRequest, UpdateRunRequest};
  use crate::runs::store::{RunStore, default_title};
  use app_commands::command::AppCommand;
  use chrono::{Local, TimeDelta, Utc};
  use helpers::{GatedReader, accept, clock_config, create, event_types, terminal_status, timetree_config};
  use maplit::btreeset;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::collections::BTreeSet;
  use std::fs;
  use std::sync::Arc;
  use std::thread;
  use tempfile::tempdir;
  use treetime::o;
  use treetime::progress::{RunWarning, RunWarningKind};
  use treetime_utils::assert_error;
  use treetime_utils::io::json::json_read_str;

  #[test]
  fn test_manager_run_records_its_lifecycle_and_outputs() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    assert_eq!(RunStatus::Created, created.status);

    let terminal = runs
      .start(&created.id, &StartRunRequest::default(), accept())
      .unwrap()
      .run();
    assert_eq!("ok", terminal_status(&terminal));

    let record = runs.get(&created.id).unwrap();
    let out_dir = runs.store().out_dir(&created.id);
    assert_eq!(
      (RunStatus::Ok, json!(out_dir), true, true, vec!["started", "terminal"],),
      (
        record.status,
        Value::Object(record.config.settings().unwrap())["output_all"].clone(),
        record
          .output_files
          .iter()
          .all(|file| out_dir.join(&file.path).is_file()),
        record.headline.clock_rate.is_some() && record.headline.r_squared.is_some(),
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
    config["output_selection"] = json!(["clock-model"]);
    config["output_all"] = json!("/elsewhere");
    config["output_clock_model"] = json!("/elsewhere/model.json");
    let created = create(&runs, AppCommand::Clock, config);
    runs
      .start(&created.id, &StartRunRequest::default(), accept())
      .unwrap()
      .run();
    let record = runs.get(&created.id).unwrap();
    let files: Vec<&str> = record
      .output_files
      .iter()
      .map(|file| file.path.to_str().unwrap())
      .collect();
    assert_eq!(
      (
        json!(["clock-model", "auspice"]),
        json!(null),
        vec!["clock.auspice.json", "clock.clock-model.json"]
      ),
      (
        Value::Object(record.config.settings().unwrap())["output_selection"].clone(),
        Value::Object(record.config.settings().unwrap())["output_clock_model"].clone(),
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
    let started = runs.start(&second.id, &StartRunRequest::default(), accept()).unwrap();
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
      runs.start(&created.id, &StartRunRequest::default(), accept()),
      format!("run `{}` has already started", created.id.as_str())
    );
  }

  #[test]
  fn test_manager_two_concurrent_runs_cancel_independently() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let first = create(&runs, AppCommand::Timetree, timetree_config());
    let second = create(&runs, AppCommand::Timetree, timetree_config());
    let first_started = runs.start(&first.id, &StartRunRequest::default(), accept()).unwrap();
    let second_started = runs.start(&second.id, &StartRunRequest::default(), accept()).unwrap();

    let cancelled = Arc::new(Mutex::new(false));
    let manager = Arc::clone(&runs);
    let first_id = first.id.clone();
    let flag = Arc::clone(&cancelled);
    runs
      .subscribe(
        &first.id,
        0,
        Box::new(move |event: &RunEvent| {
          if matches!(event.event, JobEvent::Progress { .. }) && !*flag.lock() {
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
        .filter(|event| matches!(event.event, JobEvent::Progress { .. }))
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
      let started = runs.start(&created.id, &StartRunRequest::default(), accept()).unwrap();
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
      JobEvent::Terminal {
        data: TerminalEvent::Interrupted { .. }
      }
    ));
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
      .upload_input(&created.id, "tree.nwk", &[], &mut &[b'x'; 700][..], 1024)
      .unwrap();
    assert_eq!(
      (700, runs.store().inputs_dir(&created.id).join("1/tree.nwk")),
      (uploaded.size, uploaded.path)
    );
    assert_error!(
      runs.upload_input(&created.id, "metadata.tsv", &[], &mut &[b'x'; 400][..], 1024),
      "the inputs of a run are limited to 1 KiB in total; `metadata.tsv` does not fit"
    );
    assert_eq!(
      (false, 324),
      (
        runs.store().inputs_dir(&created.id).join("2").exists(),
        runs
          .upload_input(&created.id, "metadata.tsv", &[], &mut &[b'x'; 324][..], 1024)
          .unwrap()
          .size
      )
    );
  }

  #[test]
  fn test_manager_upload_keeps_files_with_the_same_name_apart() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Ancestral, json!({}));
    let ha = runs
      .upload_input(&created.id, "aligned.fasta", &[], &mut &b">A\nACGT\n"[..], 1024)
      .unwrap();
    let na = runs
      .upload_input(&created.id, "aligned.fasta", &[], &mut &b">A\nAC\n"[..], 1024)
      .unwrap();
    let inputs = runs.store().inputs_dir(&created.id);
    assert_eq!(
      (
        (inputs.join("1/aligned.fasta"), o!(">A\nACGT\n")),
        (inputs.join("2/aligned.fasta"), o!(">A\nAC\n")),
        (o!("aligned.fasta"), o!("aligned.fasta")),
      ),
      (
        (ha.path.clone(), fs::read_to_string(&ha.path).unwrap()),
        (na.path.clone(), fs::read_to_string(&na.path).unwrap()),
        (ha.name, na.name),
      )
    );
  }

  #[test]
  fn test_manager_upload_rejects_names_with_directories() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    for name in ["../x.nwk", "a/b.nwk", ".hidden", ".."] {
      assert_error!(
        runs.upload_input(&created.id, name, &[], &mut &b"x"[..], 1024),
        format!("upload name `{name}` must be a plain file name without directories")
      );
    }
  }

  #[test]
  fn test_manager_repeated_release_keeps_a_later_upload_with_the_same_name() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Ancestral, json!({}));
    let first = runs
      .upload_input(&created.id, "aligned.fasta", &[], &mut &b">A\nACGT\n"[..], 1024)
      .unwrap();
    runs.release_input(&created.id, &first.path).unwrap();
    let second = runs
      .upload_input(&created.id, "aligned.fasta", &[], &mut &b">A\nAC\n"[..], 1024)
      .unwrap();

    runs.release_input(&created.id, &first.path).unwrap();

    assert_eq!(
      (runs.store().inputs_dir(&created.id).join("2/aligned.fasta"), true),
      (second.path.clone(), second.path.is_file())
    );
  }

  #[test]
  fn test_manager_second_upload_waits_for_the_first_and_counts_it_toward_the_limit() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let (mut first_body, started, resume) = GatedReader::new(vec![b'x'; 600]);

    let (first, second) = thread::scope(|scope| {
      let first = scope.spawn(|| runs.upload_input(&created.id, "a.nwk", &[], &mut first_body, 1000));
      started.recv().unwrap();
      let second = scope.spawn(|| runs.upload_input(&created.id, "b.nwk", &[], &mut &[b'x'; 600][..], 1000));
      resume.send(()).unwrap();
      (first.join().unwrap(), second.join().unwrap())
    });

    assert_eq!(
      (
        Some(600),
        Some(o!(
          "the inputs of a run are limited to 1000 bytes in total; `b.nwk` does not fit"
        ))
      ),
      (
        first.ok().map(|uploaded| uploaded.size),
        second.err().map(|err| err.to_string())
      )
    );
  }

  #[test]
  fn test_manager_start_during_an_upload_is_rejected_and_the_upload_completes() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let (mut body, started, resume) = GatedReader::new(vec![b'x'; 600]);
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(clock_config())),
    };

    let (start, uploaded) = thread::scope(|scope| {
      let upload = scope.spawn(|| runs.upload_input(&created.id, "a.nwk", &[], &mut body, 1000));
      started.recv().unwrap();
      let start = runs
        .start(&created.id, &request, accept())
        .err()
        .map(|err| err.to_string());
      resume.send(()).unwrap();
      (start, upload.join().unwrap())
    });

    assert_eq!(
      (
        Some(format!(
          "run `{}` cannot start while `a.nwk` is uploaded; start it again when the upload ends",
          created.id.as_str()
        )),
        true,
        RunStatus::Created
      ),
      (
        start,
        uploaded.unwrap().path.is_file(),
        runs.get(&created.id).unwrap().status
      )
    );
  }

  #[test]
  fn test_manager_upload_to_a_missing_run_keeps_no_input_lock() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();

    let missing = JobId::random();

    assert_error!(
      runs.upload_input(&missing, "a.nwk", &[], &mut &b"(A,B);"[..], 1024),
      format!("no run with id `{}`", missing.as_str())
    );
    assert_eq!(0, runs.inputs.lock().len());
  }

  #[test]
  fn test_manager_cancel_before_start_drops_the_input_lock_and_later_uploads_are_rejected() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    runs
      .upload_input(&created.id, "a.nwk", &[], &mut &b"(A,B);"[..], 1024)
      .unwrap();

    let cancelled = runs.cancel(&created.id).unwrap();
    let late = runs.upload_input(&created.id, "b.nwk", &[], &mut &b"(A,B);"[..], 1024);

    assert_eq!(
      (
        true,
        0,
        Some(format!(
          "run `{}` has already started; inputs can be uploaded only before the run starts",
          created.id.as_str()
        ))
      ),
      (
        cancelled,
        runs.inputs.lock().len(),
        late.err().map(|err| err.to_string())
      )
    );
  }

  #[test]
  fn test_manager_cancel_during_an_upload_is_rejected_and_the_upload_completes() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let (mut body, started, resume) = GatedReader::new(vec![b'x'; 600]);

    let (cancel, uploaded) = thread::scope(|scope| {
      let upload = scope.spawn(|| runs.upload_input(&created.id, "a.nwk", &[], &mut body, 1000));
      started.recv().unwrap();
      let cancel = runs.cancel(&created.id).err().map(|err| err.to_string());
      resume.send(()).unwrap();
      (cancel, upload.join().unwrap())
    });

    assert_eq!(
      (
        Some(format!(
          "run `{}` cannot cancel while `a.nwk` is uploaded; cancel it again when the upload ends",
          created.id.as_str()
        )),
        true,
        RunStatus::Created
      ),
      (
        cancel,
        uploaded.unwrap().path.is_file(),
        runs.get(&created.id).unwrap().status
      )
    );
  }

  #[test]
  fn test_manager_start_drops_the_input_lock_and_later_uploads_are_rejected() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    runs
      .upload_input(&created.id, "a.nwk", &[], &mut &b"(A,B);"[..], 1024)
      .unwrap();
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(clock_config())),
    };
    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();

    let late = runs.upload_input(&created.id, "b.nwk", &[], &mut &b"(A,B);"[..], 1024);

    assert_eq!(
      (
        "ok",
        0,
        Some(format!(
          "run `{}` has already started; inputs can be uploaded only before the run starts",
          created.id.as_str()
        ))
      ),
      (
        terminal_status(&terminal),
        runs.inputs.lock().len(),
        late.err().map(|err| err.to_string())
      )
    );
  }

  #[test]
  fn test_manager_start_removes_the_uploads_that_the_config_does_not_use() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let zika = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
    let upload = |name: &str, bytes: Vec<u8>| {
      runs
        .upload_input(&created.id, name, &[], &mut &bytes[..], 1 << 20)
        .unwrap()
        .path
    };
    let replaced = upload("tree.nwk", b"(A,B);".to_vec());
    let tree = upload("tree.nwk", fs::read(zika.join("tree.nwk")).unwrap());
    let metadata = upload("metadata.tsv", fs::read(zika.join("metadata.tsv")).unwrap());
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(json!({ "tree": tree, "metadata": metadata }))),
    };

    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();

    let inputs = runs.store().inputs_dir(&created.id);
    let remaining = fs::read_dir(&inputs)
      .unwrap()
      .map(|entry| entry.unwrap().file_name().into_string().unwrap())
      .collect::<BTreeSet<_>>();
    assert_eq!(
      ("ok", btreeset! { o!("2"), o!("3") }, false),
      (terminal_status(&terminal), remaining, replaced.exists())
    );
  }

  #[test]
  fn test_manager_resolves_output_files_inside_the_output_folder_only() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    runs
      .start(&created.id, &StartRunRequest::default(), accept())
      .unwrap()
      .run();
    fs::write(root.path().join("outside.txt"), "secret").unwrap();
    runs.file_path(&created.id, "clock.clock-model.json").unwrap();
    for escape in ["../run.json", "/etc/passwd", "../../outside.txt", ""] {
      assert_error!(
        runs.file_path(&created.id, escape),
        format!("file path `{escape}` must name a file inside the run's output folder")
      );
    }
    assert_eq!(runs.store().out_dir(&created.id), runs.out_dir(&created.id).unwrap());
  }

  #[test]
  fn test_manager_store_writes_run_json_atomically() {
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

  #[test]
  fn test_manager_start_replaces_the_command_given_at_creation() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Timetree, json!({}));
    let request = StartRunRequest {
      command: Some(AppCommand::Clock),
      config: Some(sparse(clock_config())),
    };
    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();
    let record = runs.get(&created.id).unwrap();
    assert_eq!(
      ("ok", AppCommand::Clock),
      (terminal_status(&terminal), record.config.command())
    );
  }

  #[test]
  fn test_manager_open_removes_unstarted_runs_older_than_a_day() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let recent = create(&runs, AppCommand::Clock, clock_config());
    let mut stale = create(&runs, AppCommand::Clock, clock_config());
    let mut finished = create(&runs, AppCommand::Clock, clock_config());
    runs
      .start(&finished.id, &StartRunRequest::default(), accept())
      .unwrap()
      .run();
    finished = runs.get(&finished.id).unwrap();
    let two_days_ago = Utc::now() - TimeDelta::days(2);
    stale.created_at = two_days_ago;
    finished.created_at = two_days_ago;
    runs.store().write(&stale).unwrap();
    runs.store().write(&finished).unwrap();
    drop(runs);

    let reopened = RunManager::open(root.path()).unwrap();
    let kept: BTreeSet<JobId> = reopened.list().unwrap().runs.into_iter().map(|run| run.id).collect();
    assert_eq!(
      (btreeset! { recent.id, finished.id }, false),
      (kept, reopened.store().run_dir(&stale.id).exists())
    );
  }

  #[test]
  fn test_manager_failed_run_keeps_the_warnings_raised_before_the_failure() {
    let root = tempdir().unwrap();
    let inputs = tempdir().unwrap();
    let tree = inputs.path().join("tree.nwk");
    let metadata = inputs.path().join("metadata.tsv");
    fs::write(&tree, "((A:0.1,A:0.1)X:0.1,B:0.1)root;\n").unwrap();
    fs::write(&metadata, "strain\tdate\nZ\t2000.5\n").unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({ "tree": tree, "metadata": metadata }));

    let terminal = runs
      .start(&created.id, &StartRunRequest::default(), accept())
      .unwrap()
      .run();

    let record = runs.get(&created.id).unwrap();
    assert_eq!(
      (
        "error",
        RunStatus::Error,
        vec![RunWarning {
          kind: RunWarningKind::DuplicateNodeNames,
          message: format!(
            "The tree '{}' gives the same name to more than one node: A. Nodes with the same name receive the same data from the other inputs, and augur node data keeps one entry per name.",
            tree.display()
          ),
          names: vec![o!("A")],
        }],
      ),
      (terminal_status(&terminal), record.status, record.warnings)
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::sparse;
    use crate::job::TerminalEvent;
    use crate::runs::events::RunEvent;
    use crate::runs::manager::{ConfigHook, RunManager};
    use crate::runs::record::{CreateRunRequest, RunRecord};
    use app_commands::command::AppCommand;
    use serde_json::{Value, json};
    use std::io::{self, Read};
    use std::path::Path;
    use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
    use treetime_utils::io::json::to_json_value;

    pub(super) struct GatedReader {
      bytes: Vec<u8>,
      started: Option<SyncSender<()>>,
      resume: Receiver<()>,
    }

    impl GatedReader {
      pub(super) fn new(bytes: Vec<u8>) -> (Self, Receiver<()>, SyncSender<()>) {
        let (started_sender, started) = sync_channel(1);
        let (resume, resume_receiver) = sync_channel(1);
        let reader = Self {
          bytes,
          started: Some(started_sender),
          resume: resume_receiver,
        };
        (reader, started, resume)
      }
    }

    impl Read for GatedReader {
      fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if let Some(started) = self.started.take() {
          started.send(()).map_err(io::Error::other)?;
          self.resume.recv().map_err(io::Error::other)?;
        }
        let count = buffer.len().min(self.bytes.len());
        buffer[..count].copy_from_slice(&self.bytes[..count]);
        self.bytes.drain(..count);
        Ok(count)
      }
    }

    pub(super) fn accept() -> ConfigHook {
      Box::new(|_config: &mut Value| Ok(()))
    }

    pub(super) fn create(runs: &RunManager, command: AppCommand, config: Value) -> RunRecord {
      runs
        .create(&CreateRunRequest {
          command,
          config: sparse(config),
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
        .map(|event| match to_json_value(&event).unwrap()["type"].as_str().unwrap() {
          "started" => "started",
          "progress" => "progress",
          "log" => "log",
          "iteration" => "iteration",
          "terminal" => "terminal",
          other => panic!("unknown event type {other}"),
        })
        .collect()
    }
  }
}
