use crate::confine::PathPolicy;
use app_backend::app_paths::AppPaths;
use app_backend::app_settings::store::AppSettingsStore;
use app_backend::bridge::executor::RunExecutor;
use app_backend::bridge::service::{AppService, InputPolicy};
use app_backend::examples_download::ExampleDownloads;
use app_backend::run_limits::RunLimits;
use app_backend::runs::manager::RunManager;
use app_commands::command::AppCommand;
use eyre::Report;
use serde_json::Value;
use smart_default::SmartDefault;
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub const DEFAULT_MAX_UPLOAD_SIZE: usize = 1 << 30;

const FILES_FIELD: &str = "files";

pub fn confined_service(runs_dir: &Path, executor: Arc<dyn RunExecutor>) -> Result<Arc<AppService>, Report> {
  let runs = RunManager::open(runs_dir)?;
  let policy = Arc::new(RunsFolderInputs {
    runs: Arc::clone(&runs),
  });
  Ok(Arc::new(AppService::new(runs, policy, RunLimits::Unset, executor)))
}

pub(crate) struct AppState {
  pub config: RouterConfig,
  pub settings: Option<LocalSettings>,
  pub runs: Arc<RunManager>,
  pub service: Arc<AppService>,
}

#[derive(Clone, Debug, SmartDefault)]
pub struct RouterConfig {
  #[default(DEFAULT_MAX_UPLOAD_SIZE)]
  pub max_upload_size: usize,
}

#[derive(Clone, Debug)]
pub struct LocalSettings {
  pub store: Arc<AppSettingsStore>,
  pub paths: AppPaths,
  pub runs_error: Option<String>,
  pub examples: Arc<ExampleDownloads>,
}

struct RunsFolderInputs {
  runs: Arc<RunManager>,
}

impl InputPolicy for RunsFolderInputs {
  fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report> {
    PathPolicy::new(&self.runs.input_dirs())?.confine(command, config)
  }

  fn confine_path(&self, path: &Path) -> Result<PathBuf, Report> {
    PathPolicy::new(&self.runs.input_dirs())?.confine_path(FILES_FIELD, path)
  }
}
