#[cfg(test)]
mod tests {
  use crate::runs::app_events::{AppChange, AppEventLog};
  use crate::runs::manager::{RunManager, unconfined};
  use crate::runs::record::UpdateRunRequest;
  use helpers::{Change, changes, collect, create, recorded, seqs, summarize};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::sync::Arc;
  use tempfile::tempdir;
  use treetime_utils::o;

  #[test]
  fn test_app_events_are_numbered_from_the_first_sequence_number() {
    let log = AppEventLog::new(8, 41);
    for _ in 0..3 {
      log.append(AppChange::Resync, vec![]);
    }
    let received = collect(&log, Some(41));
    assert_eq!((vec![41, 42, 43], 43), (seqs(&received), log.head()));
  }

  #[test]
  fn test_app_events_without_from_follow_only_new_events() {
    let log = AppEventLog::new(8, 1);
    log.append(AppChange::Resync, vec![]);
    let received = collect(&log, None);
    log.append(AppChange::Resync, vec![]);
    log.append(AppChange::Resync, vec![]);
    assert_eq!(vec![2, 3], seqs(&received));
  }

  #[test]
  fn test_app_events_resume_from_a_kept_event_and_then_follow() {
    let log = AppEventLog::new(8, 1);
    for _ in 0..3 {
      log.append(AppChange::Resync, vec![]);
    }
    let received = collect(&log, Some(2));
    log.append(AppChange::Resync, vec![]);
    assert_eq!(vec![2, 3, 4], seqs(&received));
  }

  #[test]
  fn test_app_events_resume_after_the_head_sends_only_new_events() {
    let log = AppEventLog::new(8, 1);
    for _ in 0..3 {
      log.append(AppChange::Resync, vec![]);
    }
    let received = collect(&log, Some(4));
    log.append(AppChange::Resync, vec![]);
    assert_eq!(vec![4], seqs(&received));
  }

  #[test]
  fn test_app_events_resume_from_an_event_older_than_the_log_starts_with_a_resync() {
    let log = AppEventLog::new(3, 1);
    for id in ["a", "b", "c", "d", "e"] {
      log.append(AppChange::Resync, vec![format!("/api/runs/{id}")]);
    }
    let received = collect(&log, Some(2));
    log.append(AppChange::Resync, vec![]);
    assert_eq!(
      vec![
        (5, json!("resync"), json!(["/api/runs", "/api/clade-in-runs"])),
        (6, json!("resync"), json!([])),
      ],
      summarize(&received)
    );
  }

  #[test]
  fn test_app_events_resume_from_an_event_newer_than_the_head_starts_with_a_resync() {
    let log = AppEventLog::new(8, 1);
    log.append(AppChange::Resync, vec![]);
    log.append(AppChange::Resync, vec![]);
    let received = collect(&log, Some(10));
    log.append(AppChange::Resync, vec![]);
    assert_eq!(vec![2, 3], seqs(&received));
    assert_eq!(json!(["/api/runs", "/api/clade-in-runs"]), summarize(&received)[0].2);
  }

  #[test]
  fn test_app_events_resync_names_the_event_before_the_first_of_an_empty_log() {
    let log = AppEventLog::new(8, 100);
    let received = collect(&log, Some(7));
    let resumed = collect(&log, Some(100));
    log.append(AppChange::Resync, vec![]);
    assert_eq!((vec![99, 100], vec![100]), (seqs(&received), seqs(&resumed)));
  }

  #[test]
  fn test_app_events_first_sequence_number_is_at_least_one() {
    let log = AppEventLog::new(8, 0);
    let received = collect(&log, Some(5));
    log.append(AppChange::Resync, vec![]);
    assert_eq!((vec![0, 1], 1), (seqs(&received), log.head()));
  }

  #[test]
  fn test_app_events_subscriber_that_stops_receives_nothing_more() {
    let log = AppEventLog::new(8, 1);
    let received = recorded();
    let sink = Arc::clone(&received);
    log.subscribe(
      None,
      Box::new(move |event| {
        sink.lock().push(event.clone());
        false
      }),
    );
    log.append(AppChange::Resync, vec![]);
    log.append(AppChange::Resync, vec![]);
    assert_eq!(vec![1], seqs(&received));
  }

  #[test]
  fn test_app_events_record_every_change_of_a_run_with_its_stale_paths() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let received = collect(runs.app_events(), None);
    let first = create(&runs, "first");
    let second = create(&runs, "second");
    runs
      .update(
        &first.id,
        UpdateRunRequest {
          title: Some("renamed".to_owned()),
          pinned: None,
        },
      )
      .unwrap();
    runs.cancel(&second.id).unwrap();
    runs.delete(&first.id).unwrap();
    runs.restore(&first.id).unwrap();
    runs.delete(&second.id).unwrap();
    runs.purge(&second.id).unwrap();

    let stale = |id: &str| vec![o!("/api/runs"), format!("/api/runs/{id}"), o!("/api/clade-in-runs")];
    let (a, b) = (first.id.as_str(), second.id.as_str());
    assert_eq!(
      vec![
        Change::new("run-created", a, Some("first"), Some("created"), stale(a)),
        Change::new("run-created", b, Some("second"), Some("created"), stale(b)),
        Change::new("run-updated", a, Some("renamed"), Some("created"), stale(a)),
        Change::new("run-updated", b, Some("second"), Some("cancelled"), stale(b)),
        Change::new("run-deleted", a, None, None, stale(a)),
        Change::new("run-restored", a, Some("renamed"), Some("created"), stale(a)),
        Change::new("run-deleted", b, None, None, stale(b)),
        Change::new("run-purged", b, None, None, vec![]),
      ],
      changes(&received)
    );
    let numbers = seqs(&received);
    assert_eq!(
      (1..=8).map(|offset| numbers[0] + offset - 1).collect::<Vec<_>>(),
      numbers,
      "consecutive changes have consecutive numbers"
    );
  }

  #[test]
  fn test_app_events_record_the_status_changes_of_a_computed_run() {
    let root = tempdir().unwrap();
    let runs = RunManager::open(root.path()).unwrap();
    let created = create(&runs, "clock");
    let received = collect(runs.app_events(), None);
    runs.start(&created.id, None, unconfined()).unwrap().run();
    let statuses = changes(&received)
      .into_iter()
      .map(|change| (change.kind, change.status))
      .collect::<Vec<_>>();
    assert_eq!(
      vec![
        (o!("run-updated"), Some(o!("running"))),
        (o!("run-updated"), Some(o!("running"))),
        (o!("run-updated"), Some(o!("ok"))),
      ],
      statuses
    );
  }

  #[test]
  fn test_app_events_of_a_restarted_manager_are_numbered_above_the_previous_ones() {
    let root = tempdir().unwrap();
    let before = RunManager::open(root.path()).unwrap();
    let received = collect(before.app_events(), None);
    for title in ["a", "b", "c"] {
      create(&before, title);
    }
    let last = *seqs(&received).last().unwrap();
    drop(before);

    let after = RunManager::open(root.path()).unwrap();
    let resumed = collect(after.app_events(), Some(last + 1));
    assert_eq!(vec![json!("resync")], summarize(&resumed).into_iter().map(|(_, kind, _)| kind).collect::<Vec<_>>());
    assert!(after.app_events().head() >= last, "{} < {last}", after.app_events().head());
  }

  mod helpers {
    use crate::command::AppCommand;
    use crate::runs::app_events::{AppEvent, AppEventLog};
    use crate::runs::manager::RunManager;
    use crate::runs::record::{CreateRunRequest, RunRecord};
    use parking_lot::Mutex;
    use serde_json::{Value, json};
    use std::path::Path;
    use std::sync::Arc;

    pub(super) type Recorded = Arc<Mutex<Vec<AppEvent>>>;

    #[derive(Debug, PartialEq, Eq)]
    pub(super) struct Change {
      pub kind: String,
      pub id: String,
      pub title: Option<String>,
      pub status: Option<String>,
      pub stale: Vec<String>,
    }

    impl Change {
      pub(super) fn new(kind: &str, id: &str, title: Option<&str>, status: Option<&str>, stale: Vec<String>) -> Self {
        Self {
          kind: kind.to_owned(),
          id: id.to_owned(),
          title: title.map(str::to_owned),
          status: status.map(str::to_owned),
          stale,
        }
      }
    }

    pub(super) fn recorded() -> Recorded {
      Arc::new(Mutex::new(vec![]))
    }

    pub(super) fn collect(log: &AppEventLog, from: Option<usize>) -> Recorded {
      let received = recorded();
      let sink = Arc::clone(&received);
      log.subscribe(
        from,
        Box::new(move |event| {
          sink.lock().push(event.clone());
          true
        }),
      );
      received
    }

    pub(super) fn create(runs: &RunManager, title: &str) -> RunRecord {
      runs.create(request(title)).unwrap()
    }

    fn request(title: &str) -> CreateRunRequest {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      CreateRunRequest {
        command: AppCommand::Clock,
        config: json!({ "tree": zika.join("tree.nwk"), "metadata": zika.join("metadata.tsv") }),
        title: Some(title.to_owned()),
        defer_start: true,
      }
    }

    pub(super) fn seqs(received: &Recorded) -> Vec<usize> {
      received.lock().iter().map(|event| event.seq).collect()
    }

    pub(super) fn summarize(received: &Recorded) -> Vec<(usize, Value, Value)> {
      received
        .lock()
        .iter()
        .map(|event| {
          let value = serde_json::to_value(event).unwrap();
          (event.seq, value["kind"].clone(), value["stale"].clone())
        })
        .collect()
    }

    pub(super) fn changes(received: &Recorded) -> Vec<Change> {
      received
        .lock()
        .iter()
        .map(|event| {
          let value = serde_json::to_value(event).unwrap();
          let run = &value["run"];
          let text = |value: &Value| value.as_str().map(str::to_owned);
          Change {
            kind: text(&value["kind"]).unwrap(),
            id: text(&run["id"]).or_else(|| text(&value["id"])).unwrap(),
            title: text(&run["title"]),
            status: text(&run["status"]),
            stale: serde_json::from_value(value["stale"].clone()).unwrap(),
          }
        })
        .collect()
    }
  }
}
