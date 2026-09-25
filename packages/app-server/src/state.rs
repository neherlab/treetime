use crate::confine::PathPolicy;
use app_commands::command::AppCommand;
use app_commands::runs::manager::{ConfigHook, RunManager};
use eyre::Report;
use serde_json::Value;
use std::path::PathBuf;
use std::sync::Arc;

pub const DEFAULT_MAX_UPLOAD_SIZE: usize = 1 << 30;

#[derive(Clone, Debug)]
pub struct ServerConfig {
  pub data_dir: PathBuf,
  pub runs_dir: PathBuf,
  pub max_upload_size: usize,
}

pub(crate) struct AppState {
  pub config: ServerConfig,
  pub runs: Arc<RunManager>,
}

impl AppState {
  pub(crate) fn new(config: ServerConfig) -> Result<Self, Report> {
    PathPolicy::new(&config.data_dir, &[])?;
    let runs = RunManager::open(&config.runs_dir)?;
    Ok(Self { config, runs })
  }

  pub(crate) fn path_policy(&self) -> Result<PathPolicy, Report> {
    PathPolicy::new(&self.config.data_dir, &self.runs.input_dirs()?)
  }

  pub(crate) fn confine_hook(&self, command: AppCommand) -> Result<ConfigHook, Report> {
    let policy = self.path_policy()?;
    Ok(Box::new(move |config: &mut Value| policy.confine(command, config)))
  }
}
