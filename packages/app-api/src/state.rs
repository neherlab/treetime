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
use std::path::PathBuf;
use std::sync::Arc;
use tokio_util::sync::CancellationToken;
use treetime_grid::MaxGridPoints;

pub const DEFAULT_MAX_UPLOAD_SIZE: usize = 1 << 30;

pub fn server_service(config: &ServerConfig, executor: Arc<dyn RunExecutor>) -> Result<Arc<AppService>, Report> {
  let runs = RunManager::open(&config.runs_dir)?;
  let policy = Arc::new(ServerInputs {
    runs: Arc::clone(&runs),
  });
  let limits = config.max_grid_points.map_or(RunLimits::Unset, RunLimits::Server);
  Ok(Arc::new(AppService::new(runs, policy, limits, executor)))
}

pub(crate) struct AppState {
  pub config: ServerConfig,
  pub runs: Arc<RunManager>,
  pub service: Arc<AppService>,
  pub openapi: Value,
  pub instance: u64,
}

#[derive(Clone, Debug)]
pub struct ServerConfig {
  pub runs_dir: PathBuf,
  pub max_upload_size: usize,
  pub max_grid_points: Option<MaxGridPoints>,
  pub shutdown: CancellationToken,
  pub settings: Option<LocalSettings>,
}

#[derive(Clone, Debug)]
pub struct LocalSettings {
  pub store: Arc<AppSettingsStore>,
  pub paths: AppPaths,
  pub runs_error: Option<String>,
  pub examples: Arc<ExampleDownloads>,
}

struct ServerInputs {
  runs: Arc<RunManager>,
}

impl InputPolicy for ServerInputs {
  fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report> {
    PathPolicy::new(&self.runs.input_dirs()?)?.confine(command, config)
  }
}
