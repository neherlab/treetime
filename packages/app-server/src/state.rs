use crate::confine::PathPolicy;
use app_commands::bridge::service::{AppService, InputPolicy};
use app_commands::command::AppCommand;
use app_commands::runs::manager::RunManager;
use eyre::Report;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;
use tokio_util::sync::CancellationToken;

pub const DEFAULT_MAX_UPLOAD_SIZE: usize = 1 << 30;

#[derive(Clone, Debug)]
pub struct ServerConfig {
  pub data_dir: PathBuf,
  pub runs_dir: PathBuf,
  pub max_upload_size: usize,
  pub shutdown: CancellationToken,
}

pub fn server_service(config: &ServerConfig) -> Result<Arc<AppService>, Report> {
  PathPolicy::new(&config.data_dir, &[])?;
  let runs = RunManager::open(&config.runs_dir)?;
  let policy = Arc::new(ServerInputs {
    data_dir: config.data_dir.clone(),
    runs: Arc::clone(&runs),
  });
  Ok(Arc::new(AppService::new(runs, config.data_dir.clone(), policy)))
}

pub(crate) struct AppState {
  pub config: ServerConfig,
  pub runs: Arc<RunManager>,
  pub service: Arc<AppService>,
  pub openapi: Value,
  pub instance: u64,
}

struct ServerInputs {
  data_dir: PathBuf,
  runs: Arc<RunManager>,
}

impl InputPolicy for ServerInputs {
  fn confine(&self, command: AppCommand, config: &mut Value) -> Result<(), Report> {
    PathPolicy::new(&self.data_dir, &self.runs.input_dirs()?)?.confine(command, config)
  }
}
