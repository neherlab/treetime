#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::{sparse, upload_bytes};
  use crate::command::config_from_settings;
  use crate::job::{JobEvent, JobId, TerminalEvent};
  use crate::runs::events::{RunEvent, read_events};
  use crate::runs::input_store::StoredInput;
  use crate::runs::manager::RunManager;
  use crate::runs::record::{CreateRunRequest, RunRecord, RunStatus, StartRunRequest, UpdateRunRequest};
  use crate::runs::store::{NewRun, RunStore};
  use app_commands::command::AppCommand;
  use chrono::{TimeDelta, Utc};
  use helpers::{accept, allowance, clock_config, create, event_types, terminal_status, timetree_config};
  use itertools::Itertools;
  use maplit::btreeset;
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::collections::BTreeSet;
  use std::fs;
  use std::path::{Path, PathBuf};
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
      (RunStatus::Ok, json!("out"), true, true, vec!["started", "terminal"],),
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
  fn test_manager_list_skips_stored_runs_of_command_line_only_commands() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let clock = create(&runs, AppCommand::Clock, clock_config());
    let tree = json!({ "tree": "t.nwk", "alignment": ["a.fasta"] });
    runs
      .store()
      .create(NewRun {
        id: JobId::random(),
        title: "homoplasy".to_owned(),
        config: config_from_settings(AppCommand::Homoplasy, &tree).unwrap(),
        inputs: vec![],
        input_sets: vec![],
        source_run_id: None,
        example: None,
      })
      .unwrap();
    assert_eq!(
      vec![clock.id],
      runs
        .list()
        .unwrap()
        .runs
        .into_iter()
        .map(|run| run.id)
        .collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_manager_start_rejects_a_stored_run_of_a_command_line_only_command() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let tree = json!({ "tree": "t.nwk", "alignment": ["a.fasta"] });
    let record = runs
      .store()
      .create(NewRun {
        id: JobId::random(),
        title: "homoplasy".to_owned(),
        config: config_from_settings(AppCommand::Homoplasy, &tree).unwrap(),
        inputs: vec![],
        input_sets: vec![],
        source_run_id: None,
        example: None,
      })
      .unwrap();
    assert_error!(
      runs.start(&record.id, &StartRunRequest::default(), accept()),
      "the command `homoplasy` is available in the command-line tool only"
    );
  }

  #[test]
  fn test_manager_create_rejects_a_command_line_only_command() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    assert_error!(
      runs.create(&CreateRunRequest {
        command: AppCommand::Mugration,
        config: sparse(json!({ "tree": "t.nwk", "metadata": "m.tsv", "attribute": "country" })),
        defer_start: true,
        source_run_id: None,
        example: None,
      }),
      "the command `mugration` is available in the command-line tool only"
    );
  }

  #[test]
  fn test_manager_start_rejects_a_switch_to_a_command_line_only_command() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, clock_config());
    let request = StartRunRequest {
      command: Some(AppCommand::Mugration),
      config: None,
      source_run_id: None,
      example: None,
    };
    assert_error!(
      runs.start(&created.id, &request, accept()).map(|_| ()),
      "the command `mugration` is available in the command-line tool only"
    );
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
    let uploaded = upload_bytes(&runs, &created.id, "tree.nwk", &[], &[b'x'; 700], 1024).unwrap();
    assert_eq!(
      (700, runs.input_store().root().join(&uploaded.sha256).join("tree.nwk")),
      (uploaded.size, uploaded.path)
    );
    assert_error!(
      upload_bytes(&runs, &created.id, "metadata.tsv", &[], &[b'x'; 400], 1024),
      "the inputs of a run are limited to 1 KiB in total; `metadata.tsv` does not fit"
    );
    assert_eq!(
      (1, 324),
      (
        runs.input_store().keys().unwrap().len(),
        upload_bytes(&runs, &created.id, "metadata.tsv", &[], &[b'y'; 324], 1024)
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
    let ha = upload_bytes(&runs, &created.id, "aligned.fasta", &[], b">A\nACGT\n", 1024).unwrap();
    let na = upload_bytes(&runs, &created.id, "aligned.fasta", &[], b">A\nAC\n", 1024).unwrap();
    assert_eq!(
      (
        (o!(">A\nACGT\n"), o!(">A\nAC\n")),
        (o!("aligned.fasta"), o!("aligned.fasta")),
        true,
      ),
      (
        (
          fs::read_to_string(&ha.path).unwrap(),
          fs::read_to_string(&na.path).unwrap()
        ),
        (ha.name, na.name),
        ha.path != na.path,
      )
    );
  }

  #[test]
  fn test_manager_run_keeps_two_uploads_with_one_name_apart_in_its_record() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Ancestral, json!({}));
    let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
    let upload = |name: &str, bytes: Vec<u8>| {
      upload_bytes(&runs, &created.id, name, &[], &bytes, 1 << 24)
        .unwrap()
        .path
    };
    let tree = upload("tree.nwk", fs::read(data.join("tree.nwk")).unwrap());
    let first = upload("aln.fasta.xz", fs::read(data.join("aln.fasta.xz")).unwrap());
    let second = upload(
      "aln.fasta.xz",
      fs::read(data.join("../../ebola/20/aln.fasta.xz")).unwrap(),
    );
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(json!({ "tree": tree, "alignment": [first, second] }))),
      source_run_id: None,
      example: None,
    };

    let started = runs.start(&created.id, &request, accept()).unwrap().record().clone();

    assert_eq!(
      (
        json!("inputs/tree.nwk"),
        json!(["inputs/1/aln.fasta.xz", "inputs/2/aln.fasta.xz"]),
      ),
      (
        Value::Object(started.config.settings().unwrap())["tree"].clone(),
        Value::Object(started.config.settings().unwrap())["alignment"].clone(),
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
        upload_bytes(&runs, &created.id, name, &[], b"x", 1024),
        format!("upload name `{name}` must be a plain file name without directories")
      );
    }
  }

  #[test]
  fn test_manager_repeated_release_keeps_a_later_upload_with_the_same_name() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Ancestral, json!({}));
    let first = upload_bytes(&runs, &created.id, "aligned.fasta", &[], b">A\nACGT\n", 1024).unwrap();
    runs.release_input(&created.id, &first.path).unwrap();
    let second = upload_bytes(&runs, &created.id, "aligned.fasta", &[], b">A\nAC\n", 1024).unwrap();

    runs.release_input(&created.id, &first.path).unwrap();

    assert_eq!((false, true), (first.path.exists(), second.path.is_file()));
  }

  #[test]
  fn test_manager_concurrent_uploads_fail_the_one_that_ends_without_room() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let mut first = runs.begin_upload(&created.id, "a.nwk", &[], 1000).unwrap();
    first.write(&[b'x'; 600]).unwrap();

    let second = upload_bytes(&runs, &created.id, "b.nwk", &[], &[b'y'; 600], 1000);
    let first = first.finish();

    assert_eq!(
      (
        Some(600),
        Some(o!(
          "the inputs of a run are limited to 1000 bytes in total; `a.nwk` does not fit"
        )),
        1,
      ),
      (
        second.ok().map(|uploaded| uploaded.size),
        first.err().map(|err| err.to_string()),
        runs.input_store().keys().unwrap().len(),
      )
    );
  }

  #[test]
  fn test_manager_upload_over_the_limit_fails_in_its_stream_and_leaves_no_file() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let mut upload = runs.begin_upload(&created.id, "a.nwk", &[], 1000).unwrap();
    upload.write(&[b'x'; 600]).unwrap();

    let refused = upload.write(&[b'x'; 600]);
    drop(upload);

    assert_error!(
      refused,
      "the inputs of a run are limited to 1000 bytes in total; `a.nwk` does not fit"
    );
    assert_eq!(
      (Vec::<PathBuf>::new(), 0),
      (
        fs::read_dir(runs.input_store().root())
          .unwrap()
          .map(|entry| entry.unwrap().path())
          .collect_vec(),
        runs.get(&created.id).unwrap().uploads.len(),
      )
    );
  }

  #[test]
  fn test_manager_upload_allowance_counts_the_inputs_of_the_run_except_the_replaced_ones() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let tree = upload_bytes(&runs, &created.id, "tree.nwk", &[], &[b'x'; 300], 1000).unwrap();
    upload_bytes(&runs, &created.id, "metadata.tsv", &[], &[b'y'; 200], 1000).unwrap();

    assert_eq!(
      ((500, 1000, o!("1000 bytes")), (800, 1000, o!("1000 bytes"))),
      (
        allowance(&runs, &created.id, &[]),
        allowance(&runs, &created.id, &[tree.path]),
      )
    );
  }

  #[test]
  fn test_manager_start_during_an_upload_is_rejected_and_the_upload_completes() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(clock_config())),
      source_run_id: None,
      example: None,
    };
    let mut upload = runs.begin_upload(&created.id, "a.nwk", &[], 1000).unwrap();
    upload.write(&[b'x'; 600]).unwrap();

    let start = runs
      .start(&created.id, &request, accept())
      .err()
      .map(|err| err.to_string());
    let uploaded = upload.finish();

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
      upload_bytes(&runs, &missing, "a.nwk", &[], b"(A,B);", 1024),
      format!("no run with id `{}`", missing.as_str())
    );
    assert_eq!(0, runs.inputs.lock().len());
  }

  #[test]
  fn test_manager_cancel_before_start_drops_the_input_lock_and_its_uploads_and_later_uploads_are_rejected() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    upload_bytes(&runs, &created.id, "a.nwk", &[], b"(A,B);", 1024).unwrap();

    let cancelled = runs.cancel(&created.id).unwrap();
    let late = upload_bytes(&runs, &created.id, "b.nwk", &[], b"(A,B);", 1024);

    assert_eq!(
      (
        true,
        0,
        Some(format!(
          "run `{}` has already started; inputs can be uploaded only before the run starts",
          created.id.as_str()
        )),
        (0, 0),
      ),
      (
        cancelled,
        runs.inputs.lock().len(),
        late.err().map(|err| err.to_string()),
        (
          runs.get(&created.id).unwrap().uploads.len(),
          runs.input_store().keys().unwrap().len()
        ),
      )
    );
  }

  #[test]
  fn test_manager_cancel_during_an_upload_is_rejected_and_the_upload_completes() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let mut upload = runs.begin_upload(&created.id, "a.nwk", &[], 1000).unwrap();
    upload.write(&[b'x'; 600]).unwrap();

    let cancel = runs.cancel(&created.id).err().map(|err| err.to_string());
    let uploaded = upload.finish();

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
    upload_bytes(&runs, &created.id, "a.nwk", &[], b"(A,B);", 1024).unwrap();
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(clock_config())),
      source_run_id: None,
      example: None,
    };
    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();

    let late = upload_bytes(&runs, &created.id, "b.nwk", &[], b"(A,B);", 1024);

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
    let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
    let upload = |name: &str, bytes: Vec<u8>| {
      upload_bytes(&runs, &created.id, name, &[], &bytes, 1 << 20)
        .unwrap()
        .path
    };
    let replaced = upload("tree.nwk", b"(A,B);".to_vec());
    let tree = upload("tree.nwk", fs::read(zika.join("tree.nwk")).unwrap());
    let metadata = upload("metadata.tsv", fs::read(zika.join("metadata.tsv")).unwrap());
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(json!({ "tree": tree, "metadata": metadata }))),
      source_run_id: None,
      example: None,
    };

    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();

    assert_eq!(
      ("ok", false, true, true, Vec::<StoredInput>::new()),
      (
        terminal_status(&terminal),
        replaced.exists(),
        tree.is_file(),
        metadata.is_file(),
        runs.get(&created.id).unwrap().uploads
      )
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
  }

  #[test]
  fn test_manager_store_writes_run_json_atomically() {
    let root = tempdir().unwrap();
    let store = RunStore::open(root.path()).unwrap();
    let record = store
      .create(NewRun {
        id: JobId::random(),
        title: "prune".to_owned(),
        config: config_from_settings(AppCommand::Prune, &json!({ "tree": "t.nwk" })).unwrap(),
        inputs: vec![],
        input_sets: vec![],
        source_run_id: None,
        example: None,
      })
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
      (vec!["out", "run.json"], "prune"),
      (
        names.iter().map(String::as_str).collect::<Vec<_>>(),
        reread.title.as_str()
      )
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
      source_run_id: None,
      example: None,
    };
    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();
    let record = runs.get(&created.id).unwrap();
    assert_eq!(
      ("ok", AppCommand::Clock),
      (terminal_status(&terminal), record.config.command())
    );
  }

  #[test]
  fn test_manager_start_names_the_run_after_the_example_given_at_start() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, AppCommand::Clock, json!({}));
    let request = StartRunRequest {
      command: None,
      config: Some(sparse(clock_config())),
      source_run_id: None,
      example: Some(o!("zika/20")),
    };
    let terminal = runs.start(&created.id, &request, accept()).unwrap().run();
    let record = runs.get(&created.id).unwrap();
    assert_eq!(
      ("ok", Some(o!("zika/20")), true),
      (
        terminal_status(&terminal),
        record.example,
        record.title.starts_with("zika/20, ")
      )
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
            "The tree '{}' gives the same name to more than one node: A. TreeTime keeps the name of the first such node and writes the others as <name>_2, <name>_3 and so on. Nodes with the same input name still receive the same data from the other inputs.",
            tree.display()
          ),
          names: vec![o!("A")],
        }],
      ),
      (terminal_status(&terminal), record.status, record.warnings)
    );
  }

  #[test]
  fn test_manager_failed_run_names_the_settings_of_its_error() {
    let root = tempdir().unwrap();
    let inputs = tempdir().unwrap();
    let tree = inputs.path().join("tree.nwk");
    let metadata = inputs.path().join("metadata.tsv");
    fs::write(&tree, "((A:0.1,B:0.1)X:0.1,C:0.1)root;\n").unwrap();
    fs::write(&metadata, "strain\tdate\nA\t2000.5\nB\t2001.5\nC\t2002.5\n").unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(
      &runs,
      AppCommand::Timetree,
      json!({ "tree": tree, "metadata": metadata }),
    );

    runs
      .start(&created.id, &StartRunRequest::default(), accept())
      .unwrap()
      .run();

    let record = runs.get(&created.id).unwrap();
    assert_eq!(
      (RunStatus::Error, Some(vec![o!("sequence_length")])),
      (record.status, record.error.map(|error| error.settings))
    );
  }

  mod helpers {
    use crate::__tests__::test_support::tests::sparse;
    use crate::job::{JobId, TerminalEvent};
    use crate::runs::events::RunEvent;
    use crate::runs::manager::{ConfigHook, RunManager};
    use crate::runs::record::{CreateRunRequest, RunRecord};
    use app_commands::command::AppCommand;
    use serde_json::{Value, json};
    use std::path::{Path, PathBuf};
    use treetime_utils::io::json::to_json_value;

    pub(super) fn allowance(runs: &RunManager, id: &JobId, replaces: &[PathBuf]) -> (usize, usize, String) {
      let allowance = runs.upload_allowance(id, replaces, 1000).unwrap();
      (allowance.remaining, allowance.limit, allowance.limit_text)
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
          source_run_id: None,
          example: None,
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
