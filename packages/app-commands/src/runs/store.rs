use crate::atomic_write::write_atomically;
use crate::command_config::CommandConfig;
use crate::job::{JobEvent, JobId, TerminalEvent};
use crate::runs::errors::not_found;
use crate::runs::events::EventLog;
use crate::runs::headline::RunHeadline;
use crate::runs::record::{RunRecord, RunStatus};
use chrono::{DateTime, Local, TimeZone, Utc};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::warn;
use std::cmp::Reverse;
use std::fmt::Display;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use treetime_schema::version_info;
use treetime_utils::error::report_to_string;
use treetime_utils::io::fs::read_file_to_string_if_exists;
use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};

const RUN_FILE: &str = "run.json";
const EVENTS_FILE: &str = "events.jsonl";
const INPUTS_DIR: &str = "inputs";
const OUT_DIR: &str = "out";
const DEFAULT_TITLE_FORMAT: &str = "Run %Y-%m-%d %H:%M";

pub fn default_title<Tz: TimeZone>(created_at: &DateTime<Tz>) -> String
where
  Tz::Offset: Display,
{
  created_at.format(DEFAULT_TITLE_FORMAT).to_string()
}

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

  pub fn inputs_dir(&self, id: &JobId) -> PathBuf {
    self.run_dir(id).join(INPUTS_DIR)
  }

  pub fn out_dir(&self, id: &JobId) -> PathBuf {
    self.run_dir(id).join(OUT_DIR)
  }

  pub fn events_path(&self, id: &JobId) -> PathBuf {
    self.run_dir(id).join(EVENTS_FILE)
  }

  pub fn create(&self, config: CommandConfig) -> Result<RunRecord, Report> {
    let id = JobId::random();
    fs::create_dir_all(&self.root)
      .wrap_err_with(|| format!("When creating the runs directory '{}'", self.root.display()))?;
    let dir = self.run_dir(&id);
    create_new_dir(&dir).wrap_err_with(|| format!("When creating the run directory '{}'", dir.display()))?;
    for sub in [INPUTS_DIR, OUT_DIR] {
      fs::create_dir_all(dir.join(sub))
        .wrap_err_with(|| format!("When creating the directory '{}'", dir.join(sub).display()))?;
    }
    let created_at = Utc::now();
    let record = RunRecord {
      id,
      title: default_title(&created_at.with_timezone(&Local)),
      config,
      status: RunStatus::Created,
      pinned: false,
      created_at,
      started_at: None,
      finished_at: None,
      duration_seconds: None,
      treetime_version: version_info().version.to_owned(),
      inputs: vec![],
      config_hash: None,
      changed_settings: vec![],
      headline: RunHeadline::default(),
      output_files: vec![],
      error: None,
    };
    self.write(&record)?;
    Ok(record)
  }

  pub fn read(&self, id: &JobId) -> Result<RunRecord, Report> {
    read_record(&self.run_dir(id).join(RUN_FILE), id)
  }

  pub fn write(&self, record: &RunRecord) -> Result<(), Report> {
    let path = self.run_dir(&record.id).join(RUN_FILE);
    let text = format!("{}\n", json_write_str(record, JsonPretty(true))?);
    write_atomically(&path, |file| Ok(file.write_all(text.as_bytes())?))
      .wrap_err_with(|| format!("When writing the run record '{}'", path.display()))
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

  pub fn input_dirs(&self) -> Result<Vec<PathBuf>, Report> {
    Ok(self.run_ids()?.iter().map(|id| self.inputs_dir(id)).collect())
  }

  pub fn remove_unstarted(&self, created_before: DateTime<Utc>) -> Result<Vec<JobId>, Report> {
    let mut removed = vec![];
    for record in self.list()? {
      if record.status != RunStatus::Created || record.created_at >= created_before {
        continue;
      }
      let dir = self.run_dir(&record.id);
      fs::remove_dir_all(&dir).wrap_err_with(|| format!("When removing the unstarted run '{}'", dir.display()))?;
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
        events.append(JobEvent::Terminal(TerminalEvent::Interrupted {
          job_id: record.id.clone(),
        }))?;
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
      let entry = entry?;
      let name = entry.file_name().to_string_lossy().into_owned();
      if !entry.file_type()?.is_dir() || !entry.path().join(RUN_FILE).is_file() {
        continue;
      }
      if let Ok(id) = JobId::parse(&name) {
        ids.push(id);
      }
    }
    Ok(ids)
  }
}

#[allow(
  clippy::create_dir,
  reason = "creating a run directory must fail when the directory exists, so two runs never share one"
)]
fn create_new_dir(dir: &Path) -> io::Result<()> {
  fs::create_dir(dir)
}

fn read_record(path: &Path, id: &JobId) -> Result<RunRecord, Report> {
  let Some(text) = read_file_to_string_if_exists(path)? else {
    return Err(not_found(format!("no run with id `{}`", id.as_str())));
  };
  json_read_str(&text).wrap_err_with(|| format!("When reading the run record '{}'", path.display()))
}
