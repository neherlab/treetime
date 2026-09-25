#[cfg(test)]
mod tests {
  use crate::job::JobEvent;
  use crate::runs::events::{EventLog, RunEvent, read_events};
  use helpers::{message, terminal};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_events_are_numbered_and_persisted() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("events.jsonl");
    let log = EventLog::open(&path).unwrap();
    log.append(message("a")).unwrap();
    log.append(message("b")).unwrap();
    log.append(terminal()).unwrap();
    let persisted: Vec<(usize, bool)> = read_events(&path, 0)
      .unwrap()
      .iter()
      .map(|event| (event.seq, matches!(event.event, JobEvent::Terminal(_))))
      .collect();
    assert_eq!(vec![(0, false), (1, false), (2, true)], persisted);
  }

  #[test]
  fn test_events_subscriber_resumes_from_an_offset_and_then_follows() {
    let dir = tempdir().unwrap();
    let log = EventLog::open(&dir.path().join("events.jsonl")).unwrap();
    for text in ["a", "b", "c"] {
      log.append(message(text)).unwrap();
    }
    let received = Arc::new(Mutex::new(vec![]));
    let sink = Arc::clone(&received);
    log.subscribe(
      1,
      Box::new(move |event: &RunEvent| {
        sink.lock().push(event.seq);
        true
      }),
    );
    log.append(message("d")).unwrap();
    log.append(terminal()).unwrap();
    assert_eq!(vec![1, 2, 3, 4], *received.lock());
  }

  #[test]
  fn test_events_subscriber_that_stops_receives_nothing_more() {
    let dir = tempdir().unwrap();
    let log = EventLog::open(&dir.path().join("events.jsonl")).unwrap();
    let received = Arc::new(Mutex::new(vec![]));
    let sink = Arc::clone(&received);
    log.subscribe(
      0,
      Box::new(move |event: &RunEvent| {
        sink.lock().push(event.seq);
        false
      }),
    );
    log.append(message("a")).unwrap();
    log.append(message("b")).unwrap();
    assert_eq!(vec![0], *received.lock());
  }

  #[test]
  fn test_events_reopened_log_continues_numbering_and_rejects_events_after_the_terminal() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("events.jsonl");
    EventLog::open(&path).unwrap().append(message("a")).unwrap();
    let log = EventLog::open(&path).unwrap();
    assert_eq!(1, log.append(terminal()).unwrap().seq);
    assert_error!(
      EventLog::open(&path).unwrap().append(message("late")),
      format!("the event log '{}' already holds a terminal event", path.display())
    );
  }

  #[test]
  fn test_events_read_from_a_missing_file_is_empty() {
    let dir = tempdir().unwrap();
    assert!(read_events(&dir.path().join("none.jsonl"), 0).unwrap().is_empty());
  }

  mod helpers {
    use crate::job::{JobEvent, JobId, TerminalEvent};
    use treetime::progress::{LogEvent, LogLevel};

    pub(super) fn message(text: &str) -> JobEvent {
      JobEvent::Log(LogEvent {
        level: LogLevel::Info,
        message: text.to_owned(),
      })
    }

    pub(super) fn terminal() -> JobEvent {
      JobEvent::Terminal(TerminalEvent::Cancelled {
        job_id: JobId::parse("run").unwrap(),
      })
    }
  }
}
