use app_commands::bridge::operations::DesktopBackend;
use app_commands::check_config::{CheckConfigRequest, CheckConfigResponse, check_config};
use app_commands::check_inputs::{CheckInputsRequest, InputFacts, check_inputs};
use app_commands::command::AppCommand;
use app_commands::job::JobId;
use app_commands::results::auspice::{AuspiceDocument, run_auspice};
use app_commands::results::clades::{CladeInRuns, CladeRequest, clade_in_runs};
use app_commands::results::compare::{RunComparison, compare_runs};
use app_commands::results::run_results::{RunResults, run_results};
use app_commands::run_config::{RunConfigRequest, RunConfigResponse, run_config};
use app_commands::runs::files::{RunFile, write_run_zip};
use app_commands::runs::manager::RunManager;
use app_commands::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use app_datasets::{DatasetCatalog, discover_datasets};
use eyre::{Report, WrapErr};
use log::{error, info};
use serde_json::Value;
use std::fs::File;
use std::io;
use std::path::Path;
use std::sync::Arc;
use std::thread;
use strum::VariantNames;
use tempfile::NamedTempFile;
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

  pub fn save_run_file(&self, id: &JobId, relative: &str, destination: &Path) -> Result<(), Report> {
    let source = self.runs.file_path(id, relative)?;
    write_atomically(destination, |file| {
      let mut reader = File::open(&source).wrap_err_with(|| format!("When opening '{}'", source.display()))?;
      io::copy(&mut reader, file)?;
      Ok(())
    })
  }

  pub fn save_run_archive(&self, id: &JobId, destination: &Path) -> Result<(), Report> {
    self.runs.get(id)?;
    let out_dir = self.runs.store().out_dir(id);
    write_atomically(destination, |file| write_run_zip(&out_dir, id.as_str(), file))
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

fn write_atomically(destination: &Path, write: impl FnOnce(&mut File) -> Result<(), Report>) -> Result<(), Report> {
  let dir = destination
    .parent()
    .filter(|dir| !dir.as_os_str().is_empty())
    .unwrap_or_else(|| Path::new("."));
  let mut file = NamedTempFile::new_in(dir).wrap_err_with(|| format!("When creating a file in '{}'", dir.display()))?;
  write(file.as_file_mut())?;
  file
    .persist(destination)
    .wrap_err_with(|| format!("When saving '{}'", destination.display()))?;
  Ok(())
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

  fn run_auspice(&self, id: JobId) -> Result<AuspiceDocument, Report> {
    run_auspice(&self.runs, &id)
  }

  fn compare_runs(&self, id: JobId, other: JobId) -> Result<RunComparison, Report> {
    compare_runs(&self.runs, &id, &other)
  }

  fn clade_in_runs(&self, request: CladeRequest) -> Result<CladeInRuns, Report> {
    clade_in_runs(&self.runs, &request)
  }
}
