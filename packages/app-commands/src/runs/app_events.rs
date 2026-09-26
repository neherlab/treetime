use crate::job::JobId;
use crate::runs::record::RunSummary;
use chrono::{DateTime, Utc};
use parking_lot::Mutex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use strum_macros::IntoStaticStr;

pub const RUNS_PATH: &str = "/api/runs";

pub const CLADE_IN_RUNS_PATH: &str = "/api/clade-in-runs";

pub fn run_stale_paths(id: &JobId) -> Vec<String> {
  vec![
    RUNS_PATH.to_owned(),
    format!("{RUNS_PATH}/{}", id.as_str()),
    CLADE_IN_RUNS_PATH.to_owned(),
  ]
}

pub fn resync_stale_paths() -> Vec<String> {
  vec![RUNS_PATH.to_owned(), CLADE_IN_RUNS_PATH.to_owned()]
}

/// Change of the app's runs, sent on the app-wide event stream.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct AppEvent {
  /// Position of the event in the app-wide event stream. Each event has the number of the previous one plus 1, and a
  /// restarted server numbers its events above those of the previous server. Subscribing from `seq + 1` resumes after
  /// the event.
  pub seq: usize,
  /// Time the event was recorded.
  #[schemars(with = "String")]
  pub time: DateTime<Utc>,
  /// REST paths whose answers the change made stale. A path covers every path below it: `/api/runs/abc` covers
  /// `/api/runs/abc/results`.
  pub stale: Vec<String>,
  #[serde(flatten)]
  pub change: AppChange,
}

/// What changed.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema, IntoStaticStr)]
#[serde(tag = "kind", rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum AppChange {
  /// A run was created.
  RunCreated {
    /// The new run.
    run: RunSummary,
  },
  /// A run changed: it started, ended, or got a new title or pinned state.
  RunUpdated {
    /// The run after the change.
    run: RunSummary,
  },
  /// A run moved to the trash.
  RunDeleted {
    /// Id of the run.
    id: JobId,
  },
  /// A run came back from the trash.
  RunRestored {
    /// The restored run.
    run: RunSummary,
  },
  /// A deleted run was removed for good.
  RunPurged {
    /// Id of the run.
    id: JobId,
  },
  /// The stream cannot continue after the requested event, because the server no longer keeps that event or the
  /// event belongs to a previous server. Every path in `stale` must be read again; the stream continues with the
  /// events after this one.
  Resync,
}

pub type AppSubscriber = Box<dyn FnMut(&AppEvent) -> bool + Send>;

pub struct AppEventLog {
  capacity: usize,
  state: Mutex<AppEventLogState>,
}

struct AppEventLogState {
  events: VecDeque<AppEvent>,
  next_seq: usize,
  subscribers: Vec<AppSubscriber>,
}

impl AppEventLog {
  pub fn new(capacity: usize, first_seq: usize) -> Self {
    Self {
      capacity,
      state: Mutex::new(AppEventLogState {
        events: VecDeque::with_capacity(capacity),
        next_seq: first_seq.max(1),
        subscribers: vec![],
      }),
    }
  }

  pub fn append(&self, change: AppChange, stale: Vec<String>) -> AppEvent {
    let mut state = self.state.lock();
    let event = AppEvent {
      seq: state.next_seq,
      time: Utc::now(),
      stale,
      change,
    };
    state.next_seq += 1;
    state.subscribers.retain_mut(|subscriber| subscriber(&event));
    state.events.push_back(event.clone());
    while state.events.len() > self.capacity {
      state.events.pop_front();
    }
    event
  }

  pub fn subscribe(&self, from: Option<usize>, mut subscriber: AppSubscriber) {
    let mut state = self.state.lock();
    if let Some(from) = from {
      let oldest = state.events.front().map_or(state.next_seq, |event| event.seq);
      if from < oldest || from > state.next_seq {
        let resync = AppEvent {
          seq: state.next_seq - 1,
          time: Utc::now(),
          stale: resync_stale_paths(),
          change: AppChange::Resync,
        };
        if !subscriber(&resync) {
          return;
        }
      } else {
        for event in state.events.iter().filter(|event| event.seq >= from) {
          if !subscriber(event) {
            return;
          }
        }
      }
    }
    state.subscribers.push(subscriber);
  }

  pub fn head(&self) -> usize {
    self.state.lock().next_seq - 1
  }
}
