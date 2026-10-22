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
use std::cell::RefCell;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::Path;
use treetime_utils::error::report_to_string;
use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
use wasm_bindgen::JsValue;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

thread_local! {
  static EVENT_SINK: RefCell<Option<Function>> = const { RefCell::new(None) };
}

#[wasm_bindgen]
pub struct Job {
  job: RunJob,
}

#[wasm_bindgen]
impl Job {
  #[wasm_bindgen(constructor)]
  pub fn new(#[wasm_bindgen(js_name = jobJson)] job_json: &str, jobs: u32) -> Result<Self, JsError> {
    open_job(job_json, jobs)
      .map(|job| Self { job })
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

  #[wasm_bindgen(unchecked_return_type = "Promise<string>")]
  pub fn run(
    &self,
    #[wasm_bindgen(js_name = onEvent, unchecked_param_type = "(event: string) => void")] on_event: Function,
  ) -> Result<Promise, JsError> {
    EVENT_SINK.with_borrow_mut(|sink| *sink = Some(on_event));
    let outcome = self.compute();
    EVENT_SINK.with_borrow_mut(|sink| *sink = None);
    json_write_str(&outcome, JsonPretty(false))
      .map(|text| Promise::resolve(&JsValue::from(text)))
      .map_err(|report| js_error(&report))
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

  fn compute(&self) -> RunOutcome {
    let cancel = CancelToken::default();
    match run_job_step(&self.job.id, &cancel, || check_inputs(&self.job)) {
      Ok(()) => compute_run(&self.job, &cancel, |event| deliver_event(&event)),
      Err(terminal) => RunOutcome {
        terminal,
        warnings: vec![],
      },
    }
  }
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

fn deliver_event(event: &JobEvent) {
  let text = match json_write_str(event, JsonPretty(false)) {
    Ok(text) => JsValue::from(text),
    Err(err) => {
      error!("When writing an event of the run worker: {err:#}");
      return;
    },
  };
  EVENT_SINK.with_borrow(|sink| match sink {
    Some(sink) => {
      if let Err(err) = sink.call1(&JsValue::UNDEFINED, &text) {
        error!("The event receiver of the run worker failed: {err:?}");
      }
    },
    None => error!("An event of the run worker has no receiver on this thread"),
  });
}

fn js_error(report: &Report) -> JsError {
  JsError::new(&report_to_string(report))
}
