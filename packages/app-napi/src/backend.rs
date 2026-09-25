use app_commands::bridge::operations::DesktopBackend;
use app_commands::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
use app_commands::check_inputs::{CheckInputsRequest, InputFacts, check_inputs};
use app_commands::command::AppCommand;
use app_commands::job::JobId;
use app_commands::results::clades::{CladeInRuns, CladeRequest, clade_in_runs};
use app_commands::results::compare::{RunComparison, compare_runs};
use app_commands::results::run_results::{RunResults, run_results};
use app_commands::run_config::{RunConfigRequest, RunConfigResponse, run_config};
use app_commands::runs::files::RunFile;
use app_commands::runs::manager::RunManager;
use app_commands::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use app_datasets::{DatasetCatalog, discover_datasets};
use eyre::{Report, WrapErr};
use log::{error, info};
use serde_json::Value;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use strum::VariantNames;
use treetime_schema::{VersionInfo, version_info};
use treetime_utils::env::env_var_optional;

const DATA_DIR_ENV: &str = "DATA_DIR";

pub struct DesktopService {
  runs: Arc<RunManager>,
}

impl DesktopService {
  pub fn open(runs_dir: &Path) -> Result<Self, Report> {
    Ok(Self {
      runs: RunManager::open(runs_dir)?,
    })
  }

  pub fn runs(&self) -> &Arc<RunManager> {
    &self.runs
  }

  fn start(&self, id: &JobId, config: Option<Value>) -> Result<RunRecord, Report> {
    let started = self.runs.start(id, config, Box::new(|_config: &mut Value| Ok(())))?;
    let record = started.record().clone();
    let run_id = id.clone();
    thread::Builder::new()
      .name(format!("run-{}", id.as_str()))
      .spawn(move || {
        let terminal = started.run();
        match serde_json::to_string(&terminal) {
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

impl DesktopBackend for DesktopService {
  fn version(&self) -> Result<VersionInfo, Report> {
    Ok(version_info())
  }

  fn datasets(&self) -> Result<DatasetCatalog, Report> {
    let data_dir = env_var_optional(DATA_DIR_ENV)?.unwrap_or_else(|| "data".to_owned());
    discover_datasets(Path::new(&data_dir), AppCommand::VARIANTS)
  }

  fn check_config(&self, request: CheckConfigRequest) -> Result<CheckConfigResponse, Report> {
    Ok(check_config(&request))
  }

  fn run_config(&self, request: RunConfigRequest) -> Result<RunConfigResponse, Report> {
    Ok(run_config(&request, Box::new(|_config: &mut Value| Ok(()))))
  }

  fn check_inputs(&self, request: CheckInputsRequest) -> Result<InputFacts, Report> {
    Ok(check_inputs(&request))
  }

  fn list_runs(&self) -> Result<RunList, Report> {
    self.runs.list()
  }

  fn create_run(&self, request: CreateRunRequest) -> Result<RunRecord, Report> {
    let defer_start = request.defer_start;
    let record = self.runs.create(request)?;
    if defer_start {
      Ok(record)
    } else {
      self.start(&record.id, None)
    }
  }

  fn get_run(&self, id: JobId) -> Result<RunRecord, Report> {
    self.runs.get(&id)
  }

  fn start_run(&self, id: JobId, request: StartRunRequest) -> Result<RunRecord, Report> {
    self.start(&id, request.config)
  }

  fn update_run(&self, id: JobId, request: UpdateRunRequest) -> Result<RunSummary, Report> {
    self.runs.update(&id, request)
  }

  fn cancel_run(&self, id: JobId) -> Result<CancelRunResponse, Report> {
    Ok(CancelRunResponse {
      cancelled: self.runs.cancel(&id)?,
    })
  }

  fn delete_run(&self, id: JobId) -> Result<(), Report> {
    self.runs.delete(&id)
  }

  fn restore_run(&self, id: JobId) -> Result<RunSummary, Report> {
    self.runs.restore(&id)
  }

  fn purge_run(&self, id: JobId) -> Result<(), Report> {
    self.runs.purge(&id)
  }

  fn run_files(&self, id: JobId) -> Result<Vec<RunFile>, Report> {
    self.runs.files(&id)
  }

  fn run_results(&self, id: JobId) -> Result<RunResults, Report> {
    run_results(&self.runs, &id)
  }

  fn compare_runs(&self, id: JobId, other: JobId) -> Result<RunComparison, Report> {
    compare_runs(&self.runs, &id, &other)
  }

  fn clade_in_runs(&self, request: CladeRequest) -> Result<CladeInRuns, Report> {
    clade_in_runs(&self.runs, &request)
  }
}
