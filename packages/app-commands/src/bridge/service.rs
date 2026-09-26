use crate::bridge::operations::Operations;
use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
use crate::check_inputs::{CheckInputsRequest, InputFacts, check_inputs};
use crate::command::AppCommand;
use crate::datasets::{DatasetCatalog, dataset_catalog};
use crate::job::JobId;
use crate::results::auspice::{AuspiceDocument, run_auspice};
use crate::results::clades::{CladeInRuns, CladeRequest, clade_in_runs};
use crate::results::compare::{RunComparison, compare_runs};
use crate::results::run_results::{RunResults, run_results};
use crate::run_config::{RunConfigRequest, RunConfigResponse, run_config};
use crate::runs::errors::invalid;
use crate::runs::files::RunFile;
use crate::runs::manager::{ConfigHook, RunManager};
use crate::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use eyre::{Report, WrapErr};
use log::{error, info};
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;
use treetime_schema::{VersionInfo, version_info};
use treetime_utils::io::json::{JsonPretty, json_write_str};

pub trait InputPolicy: Send + Sync {
  fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report>;
}

pub struct Unconfined;

impl InputPolicy for Unconfined {
  fn confine(&self, _command: AppCommand, _config: &mut Value) -> Result<(), Report> {
    Ok(())
  }
}

pub struct AppService {
  runs: Arc<RunManager>,
  data_dir: PathBuf,
  policy: Arc<dyn InputPolicy>,
}

#[allow(
  clippy::same_name_method,
  reason = "the `Operations` implementation that dispatches an `OperationRequest` delegates to these methods"
)]
impl AppService {
  pub fn new(runs: Arc<RunManager>, data_dir: PathBuf, policy: Arc<dyn InputPolicy>) -> Self {
    Self { runs, data_dir, policy }
  }

  pub fn runs(&self) -> &Arc<RunManager> {
    &self.runs
  }

  pub fn version(&self) -> Result<VersionInfo, Report> {
    Ok(version_info())
  }

  pub fn datasets(&self) -> Result<DatasetCatalog, Report> {
    dataset_catalog(&self.data_dir)
  }

  pub fn check_config(&self, request: &CheckConfigRequest) -> Result<CheckConfigResponse, Report> {
    Ok(check_config(request))
  }

  pub fn run_config(&self, request: &RunConfigRequest) -> Result<RunConfigResponse, Report> {
    let hook = self.hook(request.command);
    Ok(run_config(request, hook))
  }

  pub fn check_inputs(&self, request: CheckInputsRequest) -> Result<InputFacts, Report> {
    let mut config = Value::Object(request.config);
    self
      .policy
      .confine(request.command, &mut config)
      .map_err(|err| invalid(format!("{err:#}")))?;
    let Value::Object(config) = config else {
      return Err(invalid("a command configuration must be a mapping of settings"));
    };
    check_inputs(&CheckInputsRequest {
      command: request.command,
      config,
    })
  }

  pub fn list_runs(&self) -> Result<RunList, Report> {
    self.runs.list()
  }

  pub fn create_run(&self, request: CreateRunRequest) -> Result<RunRecord, Report> {
    let defer_start = request.defer_start;
    let record = self.runs.create(request)?;
    if defer_start {
      Ok(record)
    } else {
      self.start(&record.id, None)
    }
  }

  pub fn get_run(&self, id: &JobId) -> Result<RunRecord, Report> {
    self.runs.get(id)
  }

  pub fn start_run(&self, id: &JobId, request: StartRunRequest) -> Result<RunRecord, Report> {
    self.start(id, request.config)
  }

  pub fn update_run(&self, id: &JobId, request: UpdateRunRequest) -> Result<RunSummary, Report> {
    self.runs.update(id, request)
  }

  pub fn cancel_run(&self, id: &JobId) -> Result<CancelRunResponse, Report> {
    Ok(CancelRunResponse {
      cancelled: self.runs.cancel(id)?,
    })
  }

  pub fn delete_run(&self, id: &JobId) -> Result<(), Report> {
    self.runs.delete(id)
  }

  pub fn restore_run(&self, id: &JobId) -> Result<RunSummary, Report> {
    self.runs.restore(id)
  }

  pub fn purge_run(&self, id: &JobId) -> Result<(), Report> {
    self.runs.purge(id)
  }

  pub fn run_files(&self, id: &JobId) -> Result<Vec<RunFile>, Report> {
    self.runs.files(id)
  }

  pub fn run_results(&self, id: &JobId) -> Result<RunResults, Report> {
    run_results(&self.runs, id)
  }

  pub fn run_auspice(&self, id: &JobId) -> Result<AuspiceDocument, Report> {
    run_auspice(&self.runs, id)
  }

  pub fn compare_runs(&self, id: &JobId, other: &JobId) -> Result<RunComparison, Report> {
    compare_runs(&self.runs, id, other)
  }

  pub fn clade_in_runs(&self, request: &CladeRequest) -> Result<CladeInRuns, Report> {
    clade_in_runs(&self.runs, request)
  }

  fn hook(&self, command: AppCommand) -> ConfigHook {
    let policy = Arc::clone(&self.policy);
    Box::new(move |config: &mut Value| policy.confine(command, config))
  }

  fn start(&self, id: &JobId, config: Option<Value>) -> Result<RunRecord, Report> {
    let command = self.runs.get(id)?.command;
    let started = self.runs.start(id, config, self.hook(command))?;
    let record = started.record().clone();
    let run_id = id.clone();
    thread::Builder::new()
      .name(format!("run-{}", id.as_str()))
      .spawn(move || {
        let terminal = started.run();
        match json_write_str(&terminal, JsonPretty(false)) {
          Ok(terminal) => info!("Run {} ended: {terminal}", run_id.as_str()),
          Err(err) => error!(
            "Run {} ended; its terminal event cannot be shown: {err}",
            run_id.as_str()
          ),
        }
      })
      .wrap_err_with(|| format!("When starting a thread for run `{}`", id.as_str()))?;
    Ok(record)
  }
}

impl Operations for AppService {
  fn version(&self) -> Result<VersionInfo, Report> {
    Self::version(self)
  }

  fn datasets(&self) -> Result<DatasetCatalog, Report> {
    Self::datasets(self)
  }

  fn check_config(&self, request: CheckConfigRequest) -> Result<CheckConfigResponse, Report> {
    Self::check_config(self, &request)
  }

  fn run_config(&self, request: RunConfigRequest) -> Result<RunConfigResponse, Report> {
    Self::run_config(self, &request)
  }

  fn check_inputs(&self, request: CheckInputsRequest) -> Result<InputFacts, Report> {
    Self::check_inputs(self, request)
  }

  fn list_runs(&self) -> Result<RunList, Report> {
    Self::list_runs(self)
  }

  fn create_run(&self, request: CreateRunRequest) -> Result<RunRecord, Report> {
    Self::create_run(self, request)
  }

  fn get_run(&self, id: JobId) -> Result<RunRecord, Report> {
    Self::get_run(self, &id)
  }

  fn start_run(&self, id: JobId, request: StartRunRequest) -> Result<RunRecord, Report> {
    Self::start_run(self, &id, request)
  }

  fn update_run(&self, id: JobId, request: UpdateRunRequest) -> Result<RunSummary, Report> {
    Self::update_run(self, &id, request)
  }

  fn cancel_run(&self, id: JobId) -> Result<CancelRunResponse, Report> {
    Self::cancel_run(self, &id)
  }

  fn delete_run(&self, id: JobId) -> Result<(), Report> {
    Self::delete_run(self, &id)
  }

  fn restore_run(&self, id: JobId) -> Result<RunSummary, Report> {
    Self::restore_run(self, &id)
  }

  fn purge_run(&self, id: JobId) -> Result<(), Report> {
    Self::purge_run(self, &id)
  }

  fn run_files(&self, id: JobId) -> Result<Vec<RunFile>, Report> {
    Self::run_files(self, &id)
  }

  fn run_results(&self, id: JobId) -> Result<RunResults, Report> {
    Self::run_results(self, &id)
  }

  fn run_auspice(&self, id: JobId) -> Result<AuspiceDocument, Report> {
    Self::run_auspice(self, &id)
  }

  fn compare_runs(&self, id: JobId, other: JobId) -> Result<RunComparison, Report> {
    Self::compare_runs(self, &id, &other)
  }

  fn clade_in_runs(&self, request: CladeRequest) -> Result<CladeInRuns, Report> {
    Self::clade_in_runs(self, &request)
  }
}
