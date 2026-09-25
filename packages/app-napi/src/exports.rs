use crate::guard::{guarded, guarded_json, to_napi};
use crate::runs::{create_run, parse_id, start_run};
use crate::subscription::EventForwarder;
use app_commands::bridge::error::ErrorResponse;
use app_commands::check_config::{CheckConfigRequest, check_config};
use app_commands::check_inputs::{CheckInputsRequest, check_inputs};
use app_commands::command::AppCommand;
use app_commands::results::clades::{CladeRequest, clade_in_runs};
use app_commands::results::compare::compare_runs;
use app_commands::results::run_results::run_results;
use app_commands::run_config::{RunConfigRequest, run_config};
use app_commands::runs::manager::RunManager;
use app_commands::runs::record::UpdateRunRequest;
use app_datasets::discover_datasets;
use eyre::Report;
use napi::bindgen_prelude::{AsyncTask, Buffer};
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use napi::{Env, Status, Task};
use napi_derive::napi;
use serde_json::Value;
use std::fs;
use std::path::Path;
use std::sync::Arc;
use strum::VariantNames;
use treetime_schema::version_info;
use treetime_utils::env::env_var_optional;
use treetime_utils::make_report;

const DATA_DIR_ENV: &str = "DATA_DIR";

type EventCallback = Arc<ThreadsafeFunction<String, ()>>;

type EventSink = Box<dyn Fn(String) -> bool + Send>;

#[napi]
pub fn version() -> napi::Result<String> {
  sync_json(|| Ok(version_info()))
}

#[napi]
pub fn datasets() -> napi::Result<String> {
  sync_json(|| {
    let data_dir = env_var_optional(DATA_DIR_ENV)?.unwrap_or_else(|| "data".to_owned());
    discover_datasets(Path::new(&data_dir), AppCommand::VARIANTS)
  })
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
pub fn check_config_json(request_json: String) -> napi::Result<String> {
  sync_json(|| {
    let request: CheckConfigRequest = serde_json::from_str(&request_json)?;
    Ok(check_config(&request))
  })
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
pub fn run_config_json(request_json: String) -> napi::Result<String> {
  sync_json(|| {
    let request: RunConfigRequest = serde_json::from_str(&request_json)?;
    Ok(run_config(&request, Box::new(|_config: &mut Value| Ok(()))))
  })
}

#[napi(ts_return_type = "Promise<string>")]
pub fn check_inputs_json(request_json: String) -> AsyncTask<JsonTask> {
  JsonTask::spawn(move || {
    let request: CheckInputsRequest = serde_json::from_str(&request_json)?;
    Ok(serde_json::to_string(&check_inputs(&request))?)
  })
}

#[napi]
pub struct RunService {
  runs: Arc<RunManager>,
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
impl RunService {
  #[napi(constructor)]
  pub fn new(runs_dir: String) -> napi::Result<Self> {
    let runs = guarded(|| RunManager::open(Path::new(&runs_dir))).map_err(|err| to_napi(&err))?;
    Ok(Self { runs })
  }

  #[napi]
  pub fn list(&self) -> napi::Result<String> {
    sync_json(|| self.runs.list())
  }

  #[napi]
  pub fn get(&self, id: String) -> napi::Result<String> {
    sync_json(|| self.runs.get(&parse_id(&id)?))
  }

  #[napi]
  pub fn create(&self, request_json: String) -> napi::Result<String> {
    sync_json(|| create_run(&self.runs, &request_json))
  }

  #[napi(
    ts_args_type = "id: string, configJson: string | null",
    ts_return_type = "Promise<string>"
  )]
  pub fn start(&self, id: String, config_json: Option<String>) -> napi::Result<AsyncTask<JsonTask>> {
    let started = guarded(|| start_run(&self.runs, &id, config_json.as_deref())).map_err(|err| to_napi(&err))?;
    Ok(JsonTask::spawn(move || Ok(serde_json::to_string(&started.run())?)))
  }

  #[napi]
  pub fn cancel(&self, id: String) -> napi::Result<bool> {
    guarded(|| self.runs.cancel(&parse_id(&id)?)).map_err(|err| to_napi(&err))
  }

  #[napi(ts_args_type = "id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void")]
  pub fn subscribe(&self, id: String, from: u32, on_event: EventCallback) -> napi::Result<Subscription> {
    let send: EventSink =
      Box::new(move |json: String| on_event.call(Ok(json), ThreadsafeFunctionCallMode::NonBlocking) == Status::Ok);
    let forwarder = EventForwarder::new(send);
    guarded(|| {
      self
        .runs
        .subscribe(&parse_id(&id)?, usize::try_from(from)?, forwarder.subscriber())
    })
    .map_err(|err| to_napi(&err))?;
    Ok(Subscription { forwarder })
  }

  #[napi]
  pub fn update(&self, id: String, request_json: String) -> napi::Result<String> {
    sync_json(|| {
      let request: UpdateRunRequest = serde_json::from_str(&request_json)?;
      self.runs.update(&parse_id(&id)?, request)
    })
  }

  #[napi]
  pub fn delete(&self, id: String) -> napi::Result<()> {
    guarded(|| self.runs.delete(&parse_id(&id)?)).map_err(|err| to_napi(&err))
  }

  #[napi]
  pub fn restore(&self, id: String) -> napi::Result<String> {
    sync_json(|| self.runs.restore(&parse_id(&id)?))
  }

  #[napi]
  pub fn purge(&self, id: String) -> napi::Result<()> {
    guarded(|| self.runs.purge(&parse_id(&id)?)).map_err(|err| to_napi(&err))
  }

  #[napi]
  pub fn files(&self, id: String) -> napi::Result<String> {
    sync_json(|| self.runs.files(&parse_id(&id)?))
  }

  #[napi]
  pub fn read_file(&self, id: String, path: String) -> napi::Result<Buffer> {
    guarded(|| {
      let path = self.runs.file_path(&parse_id(&id)?, &path)?;
      Ok(Buffer::from(fs::read(&path)?))
    })
    .map_err(|err| to_napi(&err))
  }

  #[napi]
  pub fn archive(&self, id: String) -> napi::Result<Buffer> {
    guarded(|| Ok(Buffer::from(self.runs.zip(&parse_id(&id)?)?))).map_err(|err| to_napi(&err))
  }

  #[napi(ts_return_type = "Promise<string>")]
  pub fn results(&self, id: String) -> AsyncTask<JsonTask> {
    let runs = Arc::clone(&self.runs);
    JsonTask::spawn(move || Ok(serde_json::to_string(&run_results(&runs, &parse_id(&id)?)?)?))
  }

  #[napi(ts_return_type = "Promise<string>")]
  pub fn compare(&self, id: String, other: String) -> AsyncTask<JsonTask> {
    let runs = Arc::clone(&self.runs);
    JsonTask::spawn(move || {
      let comparison = compare_runs(&runs, &parse_id(&id)?, &parse_id(&other)?)?;
      Ok(serde_json::to_string(&comparison)?)
    })
  }

  #[napi(ts_return_type = "Promise<string>")]
  pub fn clade_in_runs(&self, request_json: String) -> AsyncTask<JsonTask> {
    let runs = Arc::clone(&self.runs);
    JsonTask::spawn(move || {
      let request: CladeRequest = serde_json::from_str(&request_json)?;
      Ok(serde_json::to_string(&clade_in_runs(&runs, &request)?)?)
    })
  }
}

#[napi]
pub struct Subscription {
  forwarder: EventForwarder<EventSink>,
}

#[napi]
impl Subscription {
  #[napi]
  pub fn unsubscribe(&self) {
    self.forwarder.close();
  }
}

type Job = Box<dyn FnOnce() -> Result<String, Report> + Send>;

pub struct JsonTask {
  job: Option<Job>,
}

impl JsonTask {
  fn spawn(job: impl FnOnce() -> Result<String, Report> + Send + 'static) -> AsyncTask<Self> {
    AsyncTask::new(Self {
      job: Some(Box::new(job)),
    })
  }
}

impl Task for JsonTask {
  type Output = String;
  type JsValue = String;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    let job = self.job.take();
    guarded(|| job.map_or_else(|| Err(make_report!("the task has already run")), |job| job()))
      .map_err(|err| to_napi(&err))
  }

  fn resolve(&mut self, _env: Env, output: String) -> napi::Result<String> {
    Ok(output)
  }
}

fn sync_json<T: serde::Serialize>(operation: impl FnOnce() -> Result<T, Report>) -> napi::Result<String> {
  guarded_json(operation).map_err(|err: ErrorResponse| to_napi(&err))
}
