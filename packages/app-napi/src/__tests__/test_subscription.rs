#[cfg(test)]
mod tests {
  use crate::subscription::EventForwarder;
  use app_commands::runs::events::EventLog;
  use helpers::{started, terminal};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::sync::Arc;
  use tempfile::tempdir;

  #[test]
  fn test_subscription_forwards_events_as_json() {
    let dir = tempdir().unwrap();
    let log = EventLog::open(&dir.path().join("events.jsonl")).unwrap();
    let received = Arc::new(Mutex::new(vec![]));
    let sink = Arc::clone(&received);
    let forwarder = EventForwarder::new(move |json: String| {
      sink.lock().push(json);
      true
    });
    log.subscribe(0, forwarder.subscriber());
    log.append(started("a")).unwrap();
    let events: Vec<(Value, Value, Value)> = received
      .lock()
      .iter()
      .map(|text| {
        let event: Value = serde_json::from_str(text).unwrap();
        (event["seq"].clone(), event["type"].clone(), event["data"].clone())
      })
      .collect();
    assert_eq!(
      vec![(json!(0), json!("started"), json!({ "job_id": "a", "command": "clock" }))],
      events
    );
  }

  #[test]
  fn test_subscription_unsubscribe_stops_forwarding_and_releases_the_sink() {
    let dir = tempdir().unwrap();
    let log = EventLog::open(&dir.path().join("events.jsonl")).unwrap();
    let received = Arc::new(Mutex::new(0_usize));
    let sink = Arc::clone(&received);
    let forwarder = EventForwarder::new(move |_json: String| {
      *sink.lock() += 1;
      true
    });
    log.subscribe(0, forwarder.subscriber());
    log.append(started("a")).unwrap();
    forwarder.close();
    log.append(started("b")).unwrap();
    log.append(terminal()).unwrap();
    assert_eq!((1, 1), (*received.lock(), Arc::strong_count(&received)));
  }

  #[test]
  fn test_subscription_ends_when_the_sink_refuses_an_event() {
    let dir = tempdir().unwrap();
    let log = EventLog::open(&dir.path().join("events.jsonl")).unwrap();
    let received = Arc::new(Mutex::new(0_usize));
    let sink = Arc::clone(&received);
    let forwarder = EventForwarder::new(move |_json: String| {
      *sink.lock() += 1;
      false
    });
    log.subscribe(0, forwarder.subscriber());
    log.append(started("a")).unwrap();
    log.append(started("b")).unwrap();
    assert_eq!((1, 1), (*received.lock(), Arc::strong_count(&received)));
  }

  mod helpers {
    use app_commands::command::AppCommand;
    use app_commands::job::{JobEvent, JobId, JobStarted, TerminalEvent};

    pub(super) fn started(id: &str) -> JobEvent {
      JobEvent::Started(JobStarted {
        job_id: JobId::parse(id).unwrap(),
        command: AppCommand::Clock,
      })
    }

    pub(super) fn terminal() -> JobEvent {
      JobEvent::Terminal(TerminalEvent::Cancelled {
        job_id: JobId::parse("run").unwrap(),
      })
    }
  }
}
