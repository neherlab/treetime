use app_backend::bridge::executor::RunExecutor;
use app_backend::job::{JobEvent, JobId};
use app_backend::runs::errors::{conflict, invalid};
use app_backend::runs::files::{FileChunks, run_file_destination};
use app_backend::runs::lifecycle::{Prepared, PreparedRun, RunOutcome, StartedRun};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::Mutex;
use std::collections::BTreeMap;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::mem;
use std::path::Path;
use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
use wasm_bindgen::prelude::wasm_bindgen;

const START_ORDER: &str = "start";

const CANCEL_ORDER: &str = "cancel";

#[derive(Default)]
pub struct JobExecutor {
  queue: Mutex<JobQueue>,
}

impl JobExecutor {
  pub fn take_orders(&self) -> Result<Vec<JobOrder>, Report> {
    let mut queue = self.queue.lock();
    let orders = mem::take(&mut queue.orders);
    orders
      .into_iter()
      .filter_map(|order| match order {
        QueuedOrder::Start(id) => queue.waiting.get(&id).map(start_order),
        QueuedOrder::Cancel(id) => Some(Ok(JobOrder {
          kind: CANCEL_ORDER.to_owned(),
          run_id: id.as_str().to_owned(),
          job: None,
          input_paths: vec![],
        })),
      })
      .try_collect()
  }

  pub fn read_input(&self, run_id: &str, path: &str) -> Result<FileChunks, Report> {
    self.with_waiting(run_id, |prepared| {
      let path = Path::new(path);
      if !prepared.job().inputs.iter().any(|input| input.path == path) {
        return Err(invalid(format!(
          "'{}' is not an input file of run `{run_id}`",
          path.display()
        )));
      }
      FileChunks::open(path)
    })
  }

  pub fn append_event(&self, run_id: &str, event_json: &str) -> Result<(), Report> {
    let event: JobEvent =
      json_read_str(event_json).wrap_err_with(|| format!("When reading an event of run `{run_id}`"))?;
    self.with_waiting(run_id, |prepared| prepared.append_job_event(event))
  }

  pub fn begin_output(&self, run_id: &str, relative_path: &str) -> Result<(), Report> {
    self.with_waiting(run_id, |prepared| {
      let path = run_file_destination(&prepared.job().out_dir, relative_path)?;
      if let Some(folder) = path.parent() {
        fs::create_dir_all(folder).wrap_err_with(|| format!("When creating '{}'", folder.display()))?;
      }
      File::create(&path).wrap_err_with(|| format!("When creating '{}'", path.display()))?;
      Ok(())
    })
  }

  pub fn append_output(&self, run_id: &str, relative_path: &str, chunk: &[u8]) -> Result<(), Report> {
    self.with_waiting(run_id, |prepared| {
      let path = run_file_destination(&prepared.job().out_dir, relative_path)?;
      OpenOptions::new()
        .append(true)
        .open(&path)
        .and_then(|mut file| file.write_all(chunk))
        .wrap_err_with(|| format!("When writing '{}'", path.display()))
    })
  }

  pub fn complete(&self, run_id: &str, outcome_json: &str) -> Result<bool, Report> {
    let outcome: RunOutcome =
      json_read_str(outcome_json).wrap_err_with(|| format!("When reading the outcome of run `{run_id}`"))?;
    Ok(self.take_waiting(run_id)?.is_some_and(|prepared| {
      prepared.complete(outcome);
      true
    }))
  }

  pub fn fail(&self, run_id: &str, lines: Vec<String>) -> Result<bool, Report> {
    Ok(self.take_waiting(run_id)?.is_some_and(|prepared| {
      prepared.fail(lines);
      true
    }))
  }

  fn with_waiting<T>(&self, run_id: &str, work: impl FnOnce(&PreparedRun) -> Result<T, Report>) -> Result<T, Report> {
    let id = JobId::parse(run_id)?;
    let queue = self.queue.lock();
    let prepared = queue
      .waiting
      .get(&id)
      .ok_or_else(|| conflict(format!("run `{run_id}` is not computing in a run worker")))?;
    work(prepared)
  }

  fn take_waiting(&self, run_id: &str) -> Result<Option<PreparedRun>, Report> {
    let id = JobId::parse(run_id)?;
    Ok(self.queue.lock().waiting.remove(&id))
  }
}

impl RunExecutor for JobExecutor {
  fn execute(&self, started: StartedRun) -> Result<(), Report> {
    if let Prepared::Ready(prepared) = started.prepare() {
      let id = prepared.job().id.clone();
      let mut queue = self.queue.lock();
      queue.orders.push(QueuedOrder::Start(id.clone()));
      queue.waiting.insert(id, prepared);
    }
    Ok(())
  }

  fn cancel(&self, id: &JobId) {
    let mut queue = self.queue.lock();
    if queue.waiting.contains_key(id) {
      queue.orders.push(QueuedOrder::Cancel(id.clone()));
    }
  }
}

#[wasm_bindgen(getter_with_clone)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct JobOrder {
  pub kind: String,
  #[wasm_bindgen(js_name = runId)]
  pub run_id: String,
  pub job: Option<String>,
  #[wasm_bindgen(js_name = inputPaths)]
  pub input_paths: Vec<String>,
}

#[derive(Default)]
struct JobQueue {
  orders: Vec<QueuedOrder>,
  waiting: BTreeMap<JobId, PreparedRun>,
}

enum QueuedOrder {
  Start(JobId),
  Cancel(JobId),
}

fn start_order(prepared: &PreparedRun) -> Result<JobOrder, Report> {
  let job = prepared.job();
  Ok(JobOrder {
    kind: START_ORDER.to_owned(),
    run_id: job.id.as_str().to_owned(),
    job: Some(json_write_str(job, JsonPretty(false))?),
    input_paths: job
      .inputs
      .iter()
      .map(|input| input.path.to_string_lossy().into_owned())
      .collect_vec(),
  })
}
