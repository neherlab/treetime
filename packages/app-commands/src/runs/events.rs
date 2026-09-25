use crate::job::JobEvent;
use chrono::{DateTime, Utc};
use eyre::{Report, WrapErr};
use parking_lot::Mutex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::fs::{self, File, OpenOptions};
use std::io::{ErrorKind, Write};
use std::path::{Path, PathBuf};
use treetime_utils::make_error;

/// Event of a run, as stored in the run's `events.jsonl` and sent to subscribers.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct RunEvent {
  /// Position of the event in the run's event stream, starting at 0. Subscribing from `seq + 1` resumes after it.
  pub seq: usize,
  /// Time the event was recorded.
  #[schemars(with = "String")]
  pub time: DateTime<Utc>,
  #[serde(flatten)]
  pub event: JobEvent,
}

pub fn read_events(path: &Path, from: usize) -> Result<Vec<RunEvent>, Report> {
  let text = match fs::read_to_string(path) {
    Ok(text) => text,
    Err(err) if err.kind() == ErrorKind::NotFound => return Ok(vec![]),
    Err(err) => return Err(Report::new(err).wrap_err(format!("When reading events from '{}'", path.display()))),
  };
  text
    .lines()
    .filter(|line| !line.trim().is_empty())
    .map(|line| {
      serde_json::from_str::<RunEvent>(line)
        .wrap_err_with(|| format!("When reading an event from '{}'", path.display()))
    })
    .filter(|event| event.as_ref().map_or(true, |event| event.seq >= from))
    .collect()
}

pub type Subscriber = Box<dyn FnMut(&RunEvent) -> bool + Send>;

pub struct EventLog {
  path: PathBuf,
  state: Mutex<EventLogState>,
}

struct EventLogState {
  file: File,
  events: Vec<RunEvent>,
  subscribers: Vec<Subscriber>,
  closed: bool,
}

impl EventLog {
  pub fn open(path: &Path) -> Result<Self, Report> {
    let events = read_events(path, 0)?;
    let file = OpenOptions::new()
      .create(true)
      .append(true)
      .open(path)
      .wrap_err_with(|| format!("When opening the event log '{}'", path.display()))?;
    let closed = events
      .last()
      .is_some_and(|event| matches!(event.event, JobEvent::Terminal(_)));
    Ok(Self {
      path: path.to_path_buf(),
      state: Mutex::new(EventLogState {
        file,
        events,
        subscribers: vec![],
        closed,
      }),
    })
  }

  pub fn append(&self, event: JobEvent) -> Result<RunEvent, Report> {
    let mut state = self.state.lock();
    if state.closed {
      return make_error!("the event log '{}' already holds a terminal event", self.path.display());
    }
    let is_terminal = matches!(event, JobEvent::Terminal(_));
    let event = RunEvent {
      seq: state.events.len(),
      time: Utc::now(),
      event,
    };
    let line = format!("{}\n", serde_json::to_string(&event)?);
    state
      .file
      .write_all(line.as_bytes())
      .wrap_err_with(|| format!("When writing to the event log '{}'", self.path.display()))?;
    state.subscribers.retain_mut(|subscriber| subscriber(&event));
    state.events.push(event.clone());
    if is_terminal {
      state.closed = true;
      state.subscribers.clear();
    }
    Ok(event)
  }

  pub fn subscribe(&self, from: usize, mut subscriber: Subscriber) {
    let mut state = self.state.lock();
    for event in state.events.iter().filter(|event| event.seq >= from) {
      if !subscriber(event) {
        return;
      }
    }
    if !state.closed {
      state.subscribers.push(subscriber);
    }
  }

  pub fn is_closed(&self) -> bool {
    self.state.lock().closed
  }
}
