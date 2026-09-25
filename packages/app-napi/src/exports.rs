use crate::runs::{create_run, parse_id, start_run};
use app_commands::check_config::{CheckConfigRequest, check_config};
use app_commands::check_inputs::{CheckInputsRequest, check_inputs};
use app_commands::command::AppCommand;
use app_commands::run_config::{RunConfigRequest, run_config};
use app_commands::runs::events::RunEvent;
use app_commands::runs::manager::{RunManager, StartedRun};
use app_commands::runs::record::UpdateRunRequest;
use app_datasets::discover_datasets;
use log::error;
use napi::bindgen_prelude::{AsyncTask, Buffer};
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use napi::{Status, Task};
use napi_derive::napi;
use serde::Serialize;
use serde_json::Value;
use std::fs;
use std::path::Path;
use std::sync::Arc;
use strum::VariantNames;
use treetime_schema::version_info;
use treetime_utils::env::env_var_optional;

const DATA_DIR_ENV: &str = "DATA_DIR";

#[napi]
pub fn version() -> napi::Result<String> {
  to_json(&version_info())
}

#[napi]
pub fn datasets() -> napi::Result<String> {
  let data_dir = env_var_optional(DATA_DIR_ENV)
    .map_err(|err| to_napi(&err))?
    .unwrap_or_else(|| "data".to_owned());
  let catalog = discover_datasets(Path::new(&data_dir), AppCommand::VARIANTS).map_err(|err| to_napi(&err))?;
  to_json(&catalog)
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
pub fn check_config_json(request_json: String) -> napi::Result<String> {
  let request: CheckConfigRequest = serde_json::from_str(&request_json).map_err(|err| to_napi(&err.into()))?;
  to_json(&check_config(&request))
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
pub fn run_config_json(request_json: String) -> napi::Result<String> {
  let request: RunConfigRequest = serde_json::from_str(&request_json).map_err(|err| to_napi(&err.into()))?;
  to_json(&run_config(&request, Box::new(|_config: &mut Value| Ok(()))))
}

#[napi(ts_return_type = "Promise<string>")]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
pub fn check_inputs_json(request_json: String) -> napi::Result<AsyncTask<CheckInputsTask>> {
  let request: CheckInputsRequest = serde_json::from_str(&request_json).map_err(|err| to_napi(&err.into()))?;
  Ok(AsyncTask::new(CheckInputsTask { request }))
}

pub struct CheckInputsTask {
  request: CheckInputsRequest,
}

impl Task for CheckInputsTask {
  type Output = String;
  type JsValue = String;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    to_json(&check_inputs(&self.request))
  }

  fn resolve(&mut self, _env: napi::Env, output: String) -> napi::Result<String> {
    Ok(output)
  }
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
    Ok(Self {
      runs: RunManager::open(Path::new(&runs_dir)).map_err(|err| to_napi(&err))?,
    })
  }

  #[napi]
  pub fn list(&self) -> napi::Result<String> {
    to_json(&self.runs.list().map_err(|err| to_napi(&err))?)
  }

  #[napi]
  pub fn get(&self, id: String) -> napi::Result<String> {
    to_json(
      &self
        .runs
        .get(&parse_id(&id).map_err(|err| to_napi(&err))?)
        .map_err(|err| to_napi(&err))?,
    )
  }

  #[napi]
  pub fn create(&self, request_json: String) -> napi::Result<String> {
    to_json(&create_run(&self.runs, &request_json).map_err(|err| to_napi(&err))?)
  }

  #[napi(
    ts_args_type = "id: string, configJson: string | null",
    ts_return_type = "Promise<string>"
  )]
  pub fn start(&self, id: String, config_json: Option<String>) -> napi::Result<AsyncTask<RunTask>> {
    let started = start_run(&self.runs, &id, config_json.as_deref()).map_err(|err| to_napi(&err))?;
    Ok(AsyncTask::new(RunTask { run: Some(started) }))
  }

  #[napi]
  pub fn cancel(&self, id: String) -> napi::Result<bool> {
    self
      .runs
      .cancel(&parse_id(&id).map_err(|err| to_napi(&err))?)
      .map_err(|err| to_napi(&err))
  }

  #[napi(ts_args_type = "id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void")]
  pub fn subscribe(&self, id: String, from: u32, on_event: Arc<ThreadsafeFunction<String, ()>>) -> napi::Result<()> {
    let id = parse_id(&id).map_err(|err| to_napi(&err))?;
    let subscriber = Box::new(move |event: &RunEvent| match serde_json::to_string(event) {
      Ok(json) => on_event.call(Ok(json), ThreadsafeFunctionCallMode::NonBlocking) == Status::Ok,
      Err(err) => {
        error!("When serializing a run event: {err}");
        false
      },
    });
    self
      .runs
      .subscribe(
        &id,
        usize::try_from(from).map_err(|err| to_napi(&err.into()))?,
        subscriber,
      )
      .map_err(|err| to_napi(&err))
  }

  #[napi]
  pub fn update(&self, id: String, request_json: String) -> napi::Result<String> {
    let request: UpdateRunRequest = serde_json::from_str(&request_json).map_err(|err| to_napi(&err.into()))?;
    let id = parse_id(&id).map_err(|err| to_napi(&err))?;
    to_json(&self.runs.update(&id, request).map_err(|err| to_napi(&err))?)
  }

  #[napi]
  pub fn delete(&self, id: String) -> napi::Result<()> {
    self
      .runs
      .delete(&parse_id(&id).map_err(|err| to_napi(&err))?)
      .map_err(|err| to_napi(&err))
  }

  #[napi]
  pub fn restore(&self, id: String) -> napi::Result<String> {
    let id = parse_id(&id).map_err(|err| to_napi(&err))?;
    to_json(&self.runs.restore(&id).map_err(|err| to_napi(&err))?)
  }

  #[napi]
  pub fn purge(&self, id: String) -> napi::Result<()> {
    self
      .runs
      .purge(&parse_id(&id).map_err(|err| to_napi(&err))?)
      .map_err(|err| to_napi(&err))
  }

  #[napi]
  pub fn files(&self, id: String) -> napi::Result<String> {
    let id = parse_id(&id).map_err(|err| to_napi(&err))?;
    to_json(&self.runs.files(&id).map_err(|err| to_napi(&err))?)
  }

  #[napi]
  pub fn read_file(&self, id: String, path: String) -> napi::Result<Buffer> {
    let id = parse_id(&id).map_err(|err| to_napi(&err))?;
    let path = self.runs.file_path(&id, &path).map_err(|err| to_napi(&err))?;
    let bytes = fs::read(&path).map_err(|err| to_napi(&err.into()))?;
    Ok(Buffer::from(bytes))
  }

  #[napi]
  pub fn archive(&self, id: String) -> napi::Result<Buffer> {
    let id = parse_id(&id).map_err(|err| to_napi(&err))?;
    Ok(Buffer::from(self.runs.zip(&id).map_err(|err| to_napi(&err))?))
  }
}

pub struct RunTask {
  run: Option<StartedRun>,
}

impl Task for RunTask {
  type Output = String;
  type JsValue = String;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    let run = self
      .run
      .take()
      .ok_or_else(|| napi::Error::new(Status::GenericFailure, "the run has already been computed".to_owned()))?;
    to_json(&run.run())
  }

  fn resolve(&mut self, _env: napi::Env, output: String) -> napi::Result<String> {
    Ok(output)
  }
}

fn to_json<T: Serialize>(value: &T) -> napi::Result<String> {
  serde_json::to_string(value).map_err(|err| to_napi(&err.into()))
}

fn to_napi(err: &eyre::Report) -> napi::Error {
  napi::Error::new(Status::GenericFailure, format!("{err:#}"))
}
