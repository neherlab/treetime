use crate::atomic_write::write_atomically;
use crate::job::{JobEvent, JobId, TerminalEvent};
use crate::runs::errors::{conflict, not_found};
use crate::runs::events::EventLog;
use crate::runs::headline::RunHeadline;
use crate::runs::input_store::InputStore;
use crate::runs::record::{InputSet, RunInput, RunRecord, RunStatus};
use crate::runs::run_paths::OUT_DIR;
use app_commands::command_config::CommandConfig;
use chrono::{DateTime, Utc};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::warn;
use std::cmp::Reverse;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use tempfile::TempDir;
use treetime_schema::version_info;
use treetime_utils::error::report_to_string;
use treetime_utils::io::fs::read_file_to_string_if_exists;
use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};

pub const RUN_FILE: &str = "run.json";
pub const EVENTS_FILE: &str = "events.jsonl";
const STAGING_DIR: &str = ".staging";

#[derive(Clone, Debug)]
pub struct RunStore {
  root: PathBuf,
}

impl RunStore {
  pub fn open(root: &Path) -> Result<Self, Report> {
    fs::create_dir_all(root).wrap_err_with(|| format!("When creating the runs directory '{}'", root.display()))?;
    let root = root
      .canonicalize()
      .wrap_err_with(|| format!("When resolving the runs directory '{}'", root.display()))?;
    Ok(Self { root })
  }

  pub fn run_dir(&self, id: &JobId) -> PathBuf {
    self.root.join(id.as_str())
  }

  pub fn input_store(&self) -> InputStore {
    InputStore::new(&self.root)
  }

  pub fn out_dir(&self, id: &JobId) -> PathBuf {
    self.run_dir(id).join(OUT_DIR)
  }

  pub fn events_path(&self, id: &JobId) -> PathBuf {
    self.run_dir(id).join(EVENTS_FILE)
  }

  pub fn create(&self, run: NewRun) -> Result<RunRecord, Report> {
    let NewRun {
      id,
      title,
      config,
      inputs,
      input_sets,
      source_run_id,
      example,
    } = run;
    fs::create_dir_all(&self.root)
      .wrap_err_with(|| format!("When creating the runs directory '{}'", self.root.display()))?;
    let dir = self.run_dir(&id);
    create_new_dir(&dir).wrap_err_with(|| format!("When creating the run directory '{}'", dir.display()))?;
    fs::create_dir_all(dir.join(OUT_DIR))
      .wrap_err_with(|| format!("When creating the directory '{}'", dir.join(OUT_DIR).display()))?;
    let created_at = Utc::now();
    let record = RunRecord {
      group_id: id.clone(),
      id,
      title,
      renamed: false,
      example,
      config,
      status: RunStatus::Created,
      pinned: false,
      source_run_id,
      size: 0,
      seed: None,
      created_at,
      started_at: None,
      finished_at: None,
      duration_seconds: None,
      treetime_version: version_info().version.to_owned(),
      inputs,
      input_sets,
      uploads: vec![],
      config_hash: None,
      changed_settings: vec![],
      headline: RunHeadline::default(),
      output_files: vec![],
      warnings: vec![],
      error: None,
    };
    self.write(&record)?;
    Ok(record)
  }

  pub fn read(&self, id: &JobId) -> Result<RunRecord, Report> {
    read_record(&self.run_dir(id).join(RUN_FILE), id)
  }

  pub fn write(&self, record: &RunRecord) -> Result<(), Report> {
    write_record(&self.run_dir(&record.id), record)
  }

  pub fn contains(&self, id: &JobId) -> bool {
    self.run_dir(id).exists()
  }

  pub fn staging_dir(&self) -> Result<TempDir, Report> {
    let staging = self.root.join(STAGING_DIR);
    fs::create_dir_all(&staging).wrap_err_with(|| format!("When creating the folder '{}'", staging.display()))?;
    tempfile::Builder::new()
      .tempdir_in(&staging)
      .wrap_err_with(|| format!("When creating a folder in '{}'", staging.display()))
  }

  pub fn publish(&self, staged: &Path, record: &RunRecord) -> Result<(), Report> {
    write_record(staged, record)?;
    let dir = self.run_dir(&record.id);
    if dir.exists() {
      return Err(conflict(format!("a run with id `{}` exists", record.id.as_str())));
    }
    fs::rename(staged, &dir).wrap_err_with(|| format!("When moving the run folder to '{}'", dir.display()))
  }

  pub fn remove_staged(&self) -> Result<(), Report> {
    let staging = self.root.join(STAGING_DIR);
    match fs::remove_dir_all(&staging) {
      Err(err) if err.kind() != io::ErrorKind::NotFound => {
        Err(Report::new(err).wrap_err(format!("When removing the folder '{}'", staging.display())))
      },
      _ => Ok(()),
    }
  }

  pub fn list(&self) -> Result<Vec<RunRecord>, Report> {
    let mut records = vec![];
    for id in self.run_ids()? {
      let path = self.run_dir(&id).join(RUN_FILE);
      let Some(text) = read_file_to_string_if_exists(&path)? else {
        continue;
      };
      match json_read_str::<RunRecord>(&text) {
        Ok(record) => records.push(record),
        Err(report) => warn!(
          "Skipping the run folder '{}', whose record cannot be read: {}",
          path.parent().unwrap_or(&path).display(),
          report_to_string(&report)
        ),
      }
    }
    Ok(
      records
        .into_iter()
        .sorted_by_key(|record| (Reverse(record.created_at), record.id.clone()))
        .collect(),
    )
  }

  pub fn remove(&self, id: &JobId) -> Result<(), Report> {
    let dir = self.run_dir(id);
    fs::remove_dir_all(&dir).wrap_err_with(|| format!("When removing the run folder '{}'", dir.display()))
  }

  pub fn remove_unstarted(&self, created_before: DateTime<Utc>) -> Result<Vec<JobId>, Report> {
    let mut removed = vec![];
    for record in self.list()? {
      if record.status != RunStatus::Created || record.created_at >= created_before {
        continue;
      }
      self
        .remove(&record.id)
        .wrap_err_with(|| format!("When removing the unstarted run `{}`", record.id.as_str()))?;
      removed.push(record.id);
    }
    Ok(removed)
  }

  pub fn recover_interrupted(&self) -> Result<Vec<JobId>, Report> {
    let mut recovered = vec![];
    for mut record in self.list()? {
      if record.status != RunStatus::Running {
        continue;
      }
      let events = EventLog::open(&self.events_path(&record.id))?;
      if !events.is_closed() {
        events.append(JobEvent::Terminal {
          data: TerminalEvent::Interrupted {
            job_id: record.id.clone(),
          },
        })?;
      }
      record.status = RunStatus::Interrupted;
      record.finished_at = Some(Utc::now());
      self.write(&record)?;
      recovered.push(record.id);
    }
    Ok(recovered)
  }

  fn run_ids(&self) -> Result<Vec<JobId>, Report> {
    let mut ids = vec![];
    for entry in fs::read_dir(&self.root).wrap_err_with(|| format!("When listing '{}'", self.root.display()))? {
      let entry = entry.wrap_err_with(|| format!("When listing '{}'", self.root.display()))?;
      let name = entry.file_name().to_string_lossy().into_owned();
      let file_type = entry
        .file_type()
        .wrap_err_with(|| format!("When reading the type of '{}'", entry.path().display()))?;
      if !file_type.is_dir() || !entry.path().join(RUN_FILE).is_file() {
        continue;
      }
      match JobId::parse(&name) {
        Ok(id) => ids.push(id),
        Err(error) => warn!("Skipping the folder '{name}' in the runs folder: {error}"),
      }
    }
    Ok(ids)
  }
}

pub struct NewRun {
  pub id: JobId,
  pub title: String,
  pub config: CommandConfig,
  pub inputs: Vec<RunInput>,
  pub input_sets: Vec<InputSet>,
  pub source_run_id: Option<JobId>,
  pub example: Option<String>,
}

pub fn record_text(record: &RunRecord) -> Result<String, Report> {
  Ok(format!("{}\n", json_write_str(record, JsonPretty(true))?))
}

#[expect(
  clippy::create_dir,
  reason = "creating a run directory must fail when the directory exists, so two runs never share one"
)]
fn create_new_dir(dir: &Path) -> io::Result<()> {
  fs::create_dir(dir)
}

fn write_record(dir: &Path, record: &RunRecord) -> Result<(), Report> {
  let path = dir.join(RUN_FILE);
  let text = record_text(record)?;
  write_atomically(&path, |file| {
    file.write_all(text.as_bytes()).wrap_err("When writing the text")
  })
  .wrap_err_with(|| format!("When writing the run record '{}'", path.display()))
}

fn read_record(path: &Path, id: &JobId) -> Result<RunRecord, Report> {
  let Some(text) = read_file_to_string_if_exists(path)? else {
    return Err(not_found(format!("no run with id `{}`", id.as_str())));
  };
  json_read_str(&text).wrap_err_with(|| format!("When reading the run record '{}'", path.display()))
}
