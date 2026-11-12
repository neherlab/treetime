use crate::bridge::executor::RunExecutor;
use crate::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
use crate::check_inputs::{CheckInputsRequest, InputFacts, check_inputs};
use crate::clock_preview::{ClockPreview, ClockPreviewRequest, clock_preview};
use crate::config::catalog::command_settings;
use crate::file_kinds::{FileKinds, FileKindsRequest, file_kinds};
use crate::job::JobId;
use crate::results::auspice::{AuspiceDocument, run_auspice};
use crate::results::clades::{CladeInRuns, CladeRequest, clade_in_runs};
use crate::results::compare::{RunComparison, compare_runs};
use crate::results::run_results::{RunResults, run_results};
use crate::run_config::{RunConfigRequest, RunConfigResponse, run_config};
use crate::run_limits::RunLimits;
use crate::runs::errors::invalid;
use crate::runs::files::RunFile;
use crate::runs::inputs::{InputHashes, draft_inputs};
use crate::runs::manager::{ConfigHook, RunManager};
use crate::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use crate::runs::setting_differences::{ComparedConfig, SettingDifference, compared_inputs, setting_differences};
use app_commands::command::AppCommand;
use app_commands::json_value::SparseConfig;
use eyre::Report;
use itertools::Itertools;
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use treetime_schema::{VersionInfo, version_info};
use treetime_utils::make_error;

pub struct AppService {
  runs: Arc<RunManager>,
  policy: Arc<dyn InputPolicy>,
  limits: Arc<RunLimits>,
  executor: Arc<dyn RunExecutor>,
  input_hashes: InputHashes,
}

impl AppService {
  pub fn new(
    runs: Arc<RunManager>,
    policy: Arc<dyn InputPolicy>,
    limits: RunLimits,
    executor: Arc<dyn RunExecutor>,
  ) -> Self {
    Self {
      runs,
      policy,
      limits: Arc::new(limits),
      executor,
      input_hashes: InputHashes::default(),
    }
  }

  pub fn runs(&self) -> &Arc<RunManager> {
    &self.runs
  }

  pub fn version(&self) -> Result<VersionInfo, Report> {
    Ok(version_info())
  }

  pub fn check_config(&self, request: &CheckConfigRequest) -> Result<CheckConfigResponse, Report> {
    let mut response = check_config(request);
    if let (
      Some(run_id),
      CheckConfigResponse::Valid {
        command,
        config,
        run_differences,
        ..
      },
    ) = (&request.run_id, &mut response)
    {
      *run_differences = self.differences_from_run(run_id, *command, config)?;
    }
    Ok(response)
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
    Ok(run_config(&request, &hook, self.runs.input_store()))
  }

  pub fn check_inputs(&self, request: CheckInputsRequest) -> Result<InputFacts, Report> {
    let config = self.confined(request.command, request.config)?;
    check_inputs(&CheckInputsRequest {
      command: request.command,
      config,
    })
  }

  pub fn file_kinds(&self, request: &FileKindsRequest) -> Result<FileKinds, Report> {
    let files = request
      .files
      .iter()
      .map(|path| self.policy.confine_path(path))
      .try_collect()
      .map_err(|err: Report| invalid(format!("{err:#}")))?;
    Ok(file_kinds(&FileKindsRequest { files }))
  }

  pub fn clock_preview(&self, request: ClockPreviewRequest) -> Result<ClockPreview, Report> {
    let config = self.confined(request.command, request.config)?;
    clock_preview(&ClockPreviewRequest { config, ..request })
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

  pub fn delete_run(&self, id: &JobId) -> Result<(), Report> {
    self.runs.delete(id)
  }

  pub fn run_settings(&self, id: &JobId) -> Result<SparseConfig, Report> {
    self.runs.settings(id)
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

  fn differences_from_run(
    &self,
    id: &JobId,
    command: AppCommand,
    draft: &SparseConfig,
  ) -> Result<Option<Vec<SettingDifference>>, Report> {
    let record = self.runs.get(id)?;
    if record.config.command() != command {
      return Ok(None);
    }
    let catalog = command_settings(command)?;
    let mut draft = draft.0.clone();
    self.limits.apply(command, &mut draft)?;
    let store = self.runs.input_store();
    let draft_inputs = draft_inputs(command, &draft, |path| {
      self.input_hashes.sha256(store, &self.policy.confine_path(path)?)
    })?;
    let run_inputs = compared_inputs(&record.inputs);
    Ok(Some(setting_differences(
      &catalog,
      &ComparedConfig {
        settings: &record.config.settings()?,
        inputs: &run_inputs,
      },
      &ComparedConfig {
        settings: &draft,
        inputs: &draft_inputs,
      },
    )))
  }

  fn confined(&self, command: AppCommand, config: SparseConfig) -> Result<SparseConfig, Report> {
    let mut config = Value::Object(config.0);
    self
      .policy
      .confine(command, &mut config)
      .map_err(|err| invalid(format!("{err:#}")))?;
    let Value::Object(config) = config else {
      return Err(invalid("a command configuration must be a mapping of settings"));
    };
    Ok(SparseConfig(config))
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

  fn confine_path(&self, path: &Path) -> Result<PathBuf, Report>;
}

pub struct Unconfined;

impl InputPolicy for Unconfined {
  fn confine(&self, _command: AppCommand, _config: &mut Value) -> Result<(), Report> {
    Ok(())
  }

  fn confine_path(&self, path: &Path) -> Result<PathBuf, Report> {
    Ok(path.to_path_buf())
  }
}
