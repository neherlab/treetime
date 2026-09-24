use crate::confine::PathPolicy;
use app_commands::job::JobRegistry;
use eyre::Report;
use std::path::PathBuf;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct ServerConfig {
  pub data_dir: PathBuf,
  pub out_dir: PathBuf,
}

pub(crate) struct AppState {
  pub config: ServerConfig,
  pub jobs: Arc<JobRegistry>,
  pub paths: PathPolicy,
}

impl AppState {
  pub(crate) fn new(config: ServerConfig) -> Result<Self, Report> {
    let paths = PathPolicy::new(&config.data_dir, &[])?;
    Ok(Self {
      config,
      jobs: Arc::new(JobRegistry::default()),
      paths,
    })
  }
}
