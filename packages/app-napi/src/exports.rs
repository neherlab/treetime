use crate::jobs::{PendingJob, start_job};
use app_commands::command::{CheckConfigRequest, check_config};
use app_commands::job::{JobId, JobRegistry};
use app_datasets::discover_datasets;
use log::error;
use napi::Task;
use napi::bindgen_prelude::AsyncTask;
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use napi_derive::napi;
use std::path::Path;
use std::sync::Arc;
use treetime_schema::version_info;
use treetime_utils::env::env_var_optional;

const DATA_DIR_ENV: &str = "DATA_DIR";

#[napi]
pub fn version() -> napi::Result<String> {
  serde_json::to_string(&version_info()).map_err(|err| to_napi(&err.into()))
}

#[napi]
pub fn datasets() -> napi::Result<String> {
  let data_dir = env_var_optional(DATA_DIR_ENV)
    .map_err(|err| to_napi(&err))?
    .unwrap_or_else(|| "data".to_owned());
  let datasets = discover_datasets(Path::new(&data_dir)).map_err(|err| to_napi(&err))?;
  serde_json::to_string(&datasets).map_err(|err| to_napi(&err.into()))
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
pub fn check_config_json(request_json: String) -> napi::Result<String> {
  let request: CheckConfigRequest = serde_json::from_str(&request_json).map_err(|err| to_napi(&err.into()))?;
  serde_json::to_string(&check_config(&request)).map_err(|err| to_napi(&err.into()))
}

#[napi]
pub struct CommandRunner {
  jobs: Arc<JobRegistry>,
}

#[napi]
impl CommandRunner {
  #[napi(constructor)]
  pub fn new() -> Self {
    Self {
      jobs: Arc::new(JobRegistry::default()),
    }
  }

  #[napi(
    ts_args_type = "jobId: string, command: string, configJson: string, onEvent: (err: Error | null, eventJson: string) => void",
    ts_return_type = "Promise<string>"
  )]
  #[allow(
    clippy::needless_pass_by_value,
    reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
  )]
  pub fn run(
    &self,
    job_id: String,
    command: String,
    config_json: String,
    on_event: Arc<ThreadsafeFunction<String, ()>>,
  ) -> napi::Result<AsyncTask<CommandTask>> {
    let job = start_job(&self.jobs, &job_id, &command, &config_json).map_err(|err| to_napi(&err))?;
    Ok(AsyncTask::new(CommandTask {
      job: Some(job),
      on_event,
    }))
  }

  #[napi]
  #[allow(
    clippy::needless_pass_by_value,
    reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
  )]
  pub fn cancel(&self, job_id: String) -> bool {
    JobId::parse(&job_id).is_ok_and(|job_id| self.jobs.cancel(&job_id))
  }
}

pub struct CommandTask {
  job: Option<PendingJob>,
  on_event: Arc<ThreadsafeFunction<String, ()>>,
}

impl Task for CommandTask {
  type Output = String;
  type JsValue = String;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    let job = self
      .job
      .take()
      .ok_or_else(|| napi::Error::new(napi::Status::GenericFailure, "the job has already run".to_owned()))?;
    let on_event = Arc::clone(&self.on_event);
    let terminal = job.run(move |event| match serde_json::to_string(&event) {
      Ok(json) => {
        on_event.call(Ok(json), ThreadsafeFunctionCallMode::NonBlocking);
      },
      Err(err) => error!("When serializing a job event: {err}"),
    });
    serde_json::to_string(&terminal).map_err(|err| to_napi(&err.into()))
  }

  fn resolve(&mut self, _env: napi::Env, output: String) -> napi::Result<String> {
    Ok(output)
  }
}

fn to_napi(err: &eyre::Report) -> napi::Error {
  napi::Error::new(napi::Status::GenericFailure, format!("{err:#}"))
}
