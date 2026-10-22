use crate::chunks::ChunkReader;
use app_backend::job::{CancelToken, JobEvent, run_job_step};
use app_backend::runs::errors::invalid;
use app_backend::runs::files::{FileChunks, list_run_files, resolve_run_file};
use app_backend::runs::inputs::file_sha256;
use app_backend::runs::lifecycle::{JobInput, RunJob, RunOutcome, compute_run};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use js_sys::{Function, Promise};
use log::error;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::Path;
use treetime_utils::error::report_to_string;
use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

#[cfg(not(target_feature = "atomics"))]
use std::cell::RefCell;

#[cfg(target_feature = "atomics")]
use std::thread;
#[cfg(target_feature = "atomics")]
use tokio::sync::mpsc::{Sender, channel};
#[cfg(target_feature = "atomics")]
use treetime_utils::cast::u32_to_usize;
#[cfg(target_feature = "atomics")]
use treetime_utils::init::thread_pool::init_thread_pool;
#[cfg(target_feature = "atomics")]
use wasm_bindgen_futures::future_to_promise;

#[cfg(target_feature = "atomics")]
const RUN_THREAD_NAME: &str = "treetime-run";

#[cfg(target_feature = "atomics")]
const RUN_MESSAGE_CAPACITY: usize = 256;

#[cfg(not(target_feature = "atomics"))]
thread_local! {
  static EVENT_SINK: RefCell<Option<Function>> = const { RefCell::new(None) };
}

#[wasm_bindgen]
pub struct Job {
  job: RunJob,
  #[cfg(target_feature = "atomics")]
  jobs: usize,
}

#[wasm_bindgen]
impl Job {
  #[wasm_bindgen(constructor)]
  pub fn new(#[wasm_bindgen(js_name = jobJson)] job_json: &str, jobs: u32) -> Result<Self, JsError> {
    open_job(job_json, jobs)
      .map(|job| Self {
        job,
        #[cfg(target_feature = "atomics")]
        jobs: u32_to_usize(jobs),
      })
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = beginInput)]
  pub fn begin_input(&self, path: &str) -> Result<(), JsError> {
    self
      .input(path)
      .and_then(|input| {
        File::create(&input.path).wrap_err_with(|| format!("When creating input '{}'", input.path.display()))
      })
      .map(drop)
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = appendInput)]
  pub fn append_input(&self, path: &str, chunk: &[u8]) -> Result<(), JsError> {
    self
      .input(path)
      .and_then(|input| {
        OpenOptions::new()
          .append(true)
          .open(&input.path)
          .and_then(|mut file| file.write_all(chunk))
          .wrap_err_with(|| format!("When writing input '{}'", input.path.display()))
      })
      .map_err(|report| js_error(&report))
  }

  #[cfg(not(target_feature = "atomics"))]
  #[wasm_bindgen(unchecked_return_type = "Promise<string>")]
  pub fn run(
    &self,
    #[wasm_bindgen(js_name = onEvent, unchecked_param_type = "(event: string) => void")] on_event: Function,
  ) -> Result<Promise, JsError> {
    EVENT_SINK.with_borrow_mut(|sink| *sink = Some(on_event));
    let outcome = compute(
      &self.job,
      || Ok(()),
      |event| {
        EVENT_SINK.with_borrow(|sink| match sink {
          Some(sink) => deliver_event(sink, &event),
          None => error!("An event of the run worker has no receiver on this thread"),
        });
      },
    );
    EVENT_SINK.with_borrow_mut(|sink| *sink = None);
    outcome_text(&outcome)
      .map(|text| Promise::resolve(&JsValue::from(text)))
      .map_err(|report| js_error(&report))
  }

  #[cfg(target_feature = "atomics")]
  #[wasm_bindgen(unchecked_return_type = "Promise<string>")]
  pub fn run(
    &self,
    #[wasm_bindgen(js_name = onEvent, unchecked_param_type = "(event: string) => void")] on_event: Function,
  ) -> Result<Promise, JsError> {
    let (sender, mut receiver) = channel(RUN_MESSAGE_CAPACITY);
    let job = self.job.clone();
    let jobs = self.jobs;
    thread::Builder::new()
      .name(RUN_THREAD_NAME.to_owned())
      .spawn(move || compute_on_run_thread(&job, jobs, &sender))
      .wrap_err("When starting the run thread")
      .map_err(|report| js_error(&report))?;
    Ok(future_to_promise(async move {
      while let Some(message) = receiver.recv().await {
        match message {
          RunMessage::Event(event) => deliver_event(&on_event, &event),
          RunMessage::Done(outcome) => {
            return outcome_text(&outcome)
              .map(JsValue::from)
              .map_err(|report| js_error(&report).into());
          },
        }
      }
      Err(JsError::new("the run thread ended without an outcome").into())
    }))
  }

  #[wasm_bindgen(js_name = outputFiles)]
  pub fn output_files(&self) -> Result<Vec<String>, JsError> {
    list_run_files(&self.job.out_dir, &[])
      .map(|files| files.into_iter().map(|file| file.path).collect_vec())
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = readOutput)]
  pub fn read_output(&self, path: &str) -> Result<ChunkReader, JsError> {
    resolve_run_file(&self.job.out_dir, path)
      .and_then(|path| FileChunks::open(&path))
      .map(ChunkReader::from)
      .map_err(|report| js_error(&report))
  }

  fn input(&self, path: &str) -> Result<&JobInput, Report> {
    let path = Path::new(path);
    self.job.inputs.iter().find(|input| input.path == path).ok_or_else(|| {
      invalid(format!(
        "'{}' is not an input file of run `{}`",
        path.display(),
        self.job.id.as_str()
      ))
    })
  }
}

#[cfg(target_feature = "atomics")]
enum RunMessage {
  Event(JobEvent),
  Done(RunOutcome),
}

fn open_job(job_json: &str, jobs: u32) -> Result<RunJob, Report> {
  if jobs == 0 {
    return Err(invalid("a run worker needs at least 1 thread, got 0"));
  }
  let job: RunJob = json_read_str(job_json).wrap_err("When reading the job of a run worker")?;
  let folders = job
    .inputs
    .iter()
    .filter_map(|input| input.path.parent())
    .chain([job.out_dir.as_path()])
    .unique();
  for folder in folders {
    fs::create_dir_all(folder).wrap_err_with(|| format!("When creating '{}'", folder.display()))?;
  }
  Ok(job)
}

#[cfg(target_feature = "atomics")]
fn compute_on_run_thread(job: &RunJob, jobs: usize, sender: &Sender<RunMessage>) {
  let outcome = compute(
    job,
    || init_thread_pool(jobs).wrap_err_with(|| format!("When starting {jobs} threads for the run")),
    |event| send_run_message(sender, RunMessage::Event(event)),
  );
  send_run_message(sender, RunMessage::Done(outcome));
}

#[cfg(target_feature = "atomics")]
fn send_run_message(sender: &Sender<RunMessage>, message: RunMessage) {
  if sender.blocking_send(message).is_err() {
    error!("The run thread lost the receiver of its events");
  }
}

fn compute(
  job: &RunJob,
  start: impl FnOnce() -> Result<(), Report>,
  emit: impl Fn(JobEvent) + Send + Sync,
) -> RunOutcome {
  let cancel = CancelToken::default();
  match run_job_step(&job.id, &cancel, || {
    start()?;
    check_inputs(job)
  }) {
    Ok(()) => compute_run(job, &cancel, emit),
    Err(terminal) => RunOutcome {
      terminal,
      warnings: vec![],
    },
  }
}

fn check_inputs(job: &RunJob) -> Result<(), Report> {
  for input in &job.inputs {
    let (size, sha256) = file_sha256(&input.path)?;
    if size != input.size || sha256 != input.sha256 {
      return Err(invalid(format!(
        "input '{}' differs from the file the run prepared: expected {} bytes with SHA-256 {}, received {size} bytes with SHA-256 {sha256}",
        input.path.display(),
        input.size,
        input.sha256
      )));
    }
  }
  Ok(())
}

fn deliver_event(sink: &Function, event: &JobEvent) {
  let text = match json_write_str(event, JsonPretty(false)) {
    Ok(text) => JsValue::from(text),
    Err(err) => {
      error!("When writing an event of the run worker: {err:#}");
      return;
    },
  };
  if let Err(err) = sink.call1(&JsValue::UNDEFINED, &text) {
    error!("The event receiver of the run worker failed: {err:?}");
  }
}

fn outcome_text(outcome: &RunOutcome) -> Result<String, Report> {
  json_write_str(outcome, JsonPretty(false))
}

fn js_error(report: &Report) -> JsError {
  JsError::new(&report_to_string(report))
}
