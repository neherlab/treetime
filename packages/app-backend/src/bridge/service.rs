use crate::bridge::executor::RunExecutor;
use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
use crate::check_inputs::{CheckInputsRequest, InputFacts, check_inputs};
use crate::datasets::{DatasetCatalog, dataset_catalog};
use crate::job::JobId;
use crate::results::auspice::{AuspiceDocument, run_auspice};
use crate::results::clades::{CladeInRuns, CladeRequest, clade_in_runs};
use crate::results::compare::{RunComparison, compare_runs};
use crate::results::run_results::{RunResults, run_results};
use crate::run_config::{RunConfigRequest, RunConfigResponse, run_config};
use crate::run_limits::RunLimits;
use crate::runs::errors::invalid;
use crate::runs::files::RunFile;
use crate::runs::manager::{ConfigHook, RunManager};
use crate::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use app_commands::command::AppCommand;
use app_commands::json_value::SparseConfig;
use eyre::Report;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;
use treetime_schema::{VersionInfo, version_info};
use treetime_utils::make_error;

pub struct AppService {
  runs: Arc<RunManager>,
  examples_dir: PathBuf,
  policy: Arc<dyn InputPolicy>,
  limits: Arc<RunLimits>,
  executor: Arc<dyn RunExecutor>,
}

impl AppService {
  pub fn new(
    runs: Arc<RunManager>,
    examples_dir: PathBuf,
    policy: Arc<dyn InputPolicy>,
    limits: RunLimits,
    executor: Arc<dyn RunExecutor>,
  ) -> Self {
    Self {
      runs,
      examples_dir,
      policy,
      limits: Arc::new(limits),
      executor,
    }
  }

  pub fn runs(&self) -> &Arc<RunManager> {
    &self.runs
  }

  pub fn version(&self) -> Result<VersionInfo, Report> {
    Ok(version_info())
  }

  pub fn datasets(&self) -> Result<DatasetCatalog, Report> {
    dataset_catalog(&self.examples_dir)
  }

  pub fn check_config(&self, request: &CheckConfigRequest) -> Result<CheckConfigResponse, Report> {
    Ok(check_config(request))
  }

  pub fn run_config(&self, request: &RunConfigRequest) -> Result<RunConfigResponse, Report> {
    let mut config = request.config.clone();
    if let Err(report) = self.limits.apply(request.command, &mut config) {
      return Ok(RunConfigResponse::invalid(&report));
    }
    let request = RunConfigRequest {
      command: request.command,
      config,
    };
    let hook = self.hook(request.command);
    Ok(run_config(&request, hook))
  }

  pub fn check_inputs(&self, request: CheckInputsRequest) -> Result<InputFacts, Report> {
    let mut config = Value::Object(request.config.0);
    self
      .policy
      .confine(request.command, &mut config)
      .map_err(|err| invalid(format!("{err:#}")))?;
    let Value::Object(config) = config else {
      return Err(invalid("a command configuration must be a mapping of settings"));
    };
    check_inputs(&CheckInputsRequest {
      command: request.command,
      config: SparseConfig(config),
    })
  }

  pub fn list_runs(&self) -> Result<RunList, Report> {
    self.runs.list()
  }

  pub fn create_run(&self, request: &CreateRunRequest) -> Result<RunRecord, Report> {
    let defer_start = request.defer_start;
    let record = self.runs.create(request)?;
    if defer_start {
      Ok(record)
    } else {
      self.start(&record.id, &StartRunRequest::default())
    }
  }

  pub fn get_run(&self, id: &JobId) -> Result<RunRecord, Report> {
    self.runs.get(id)
  }

  pub fn start_run(&self, id: &JobId, request: &StartRunRequest) -> Result<RunRecord, Report> {
    self.start(id, request)
  }

  pub fn update_run(&self, id: &JobId, request: UpdateRunRequest) -> Result<RunSummary, Report> {
    self.runs.update(id, request)
  }

  pub fn cancel_run(&self, id: &JobId) -> Result<CancelRunResponse, Report> {
    let cancelled = self.runs.cancel(id)?;
    self.executor.cancel(id);
    Ok(CancelRunResponse { cancelled })
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
    let limits = Arc::clone(&self.limits);
    Box::new(move |config: &mut Value| {
      policy.confine(command, config)?;
      let Value::Object(settings) = config else {
        return make_error!("a command configuration must be a mapping of settings");
      };
      limits.apply(command, settings)
    })
  }

  fn start(&self, id: &JobId, request: &StartRunRequest) -> Result<RunRecord, Report> {
    let command = match request.command {
      Some(command) => command,
      None => self.runs.get(id)?.config.command(),
    };
    let started = self.runs.start(id, request, self.hook(command))?;
    let record = started.record().clone();
    self.executor.execute(started)?;
    Ok(record)
  }
}

pub trait InputPolicy: Send + Sync {
  fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report>;
}

pub struct Unconfined;

impl InputPolicy for Unconfined {
  fn confine(&self, _command: AppCommand, _config: &mut Value) -> Result<(), Report> {
    Ok(())
  }
}
