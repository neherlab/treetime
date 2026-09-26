use app_commands::bridge::service::{AppService, Unconfined};
use app_commands::job::JobId;
use app_commands::runs::files::write_run_zip;
use app_commands::runs::manager::RunManager;
use eyre::{Report, WrapErr};
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use tempfile::NamedTempFile;
use treetime_utils::env::env_var_optional;

const DATA_DIR_ENV: &str = "DATA_DIR";

const DEFAULT_DATA_DIR: &str = "data";

pub struct DesktopService {
  app: AppService,
}

impl DesktopService {
  pub fn open(runs_dir: &Path) -> Result<Self, Report> {
    let data_dir = env_var_optional(DATA_DIR_ENV)?.unwrap_or_else(|| DEFAULT_DATA_DIR.to_owned());
    Ok(Self {
      app: AppService::new(
        RunManager::open(runs_dir)?,
        PathBuf::from(data_dir),
        Arc::new(Unconfined),
      ),
    })
  }

  pub fn app(&self) -> &AppService {
    &self.app
  }

  pub fn save_run_file(&self, id: &JobId, relative: &str, destination: &Path) -> Result<(), Report> {
    let source = self.app.runs().file_path(id, relative)?;
    write_atomically(destination, |file| {
      let mut reader = File::open(&source).wrap_err_with(|| format!("When opening '{}'", source.display()))?;
      io::copy(&mut reader, file)?;
      Ok(())
    })
  }

  pub fn save_run_archive(&self, id: &JobId, destination: &Path) -> Result<(), Report> {
    let runs = self.app.runs();
    runs.get(id)?;
    let out_dir = runs.store().out_dir(id);
    write_atomically(destination, |file| write_run_zip(&out_dir, id.as_str(), file))
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
