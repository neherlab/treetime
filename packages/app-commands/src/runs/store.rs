use crate::command::AppCommand;
use crate::job::{JobEvent, JobId, TerminalEvent};
use crate::runs::errors::{conflict, invalid, not_found};
use crate::runs::events::EventLog;
use crate::runs::record::{RunRecord, RunStatus};
use chrono::Utc;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use serde_json::Value;
use std::cmp::Reverse;
use std::collections::BTreeMap;
use std::fs;
use std::io::{self, ErrorKind, Write};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;
use treetime_schema::version_info;
use treetime_utils::io::json::{JsonPretty, json_write_str};

const RUN_FILE: &str = "run.json";
const EVENTS_FILE: &str = "events.jsonl";
const INPUTS_DIR: &str = "inputs";
const OUT_DIR: &str = "out";
const TRASH_DIR: &str = ".trash";

#[derive(Clone, Debug)]
pub struct RunStore {
  root: PathBuf,
}

impl RunStore {
  pub fn open(root: &Path) -> Result<Self, Report> {
    fs::create_dir_all(root.join(TRASH_DIR))
      .wrap_err_with(|| format!("When creating the runs directory '{}'", root.display()))?;
    let root = root
      .canonicalize()
      .wrap_err_with(|| format!("When resolving the runs directory '{}'", root.display()))?;
    Ok(Self { root })
  }

  pub fn root(&self) -> &Path {
    &self.root
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

  pub fn create(&self, command: AppCommand, config: Value, title: Option<String>) -> Result<RunRecord, Report> {
    let Value::Object(config) = config else {
      return Err(invalid("a command configuration must be a mapping of settings"));
    };
    let id = JobId::random();
    let dir = self.run_dir(&id);
    create_new_dir(&dir).wrap_err_with(|| format!("When creating the run directory '{}'", dir.display()))?;
    for sub in [INPUTS_DIR, OUT_DIR] {
      fs::create_dir_all(dir.join(sub))
        .wrap_err_with(|| format!("When creating the directory '{}'", dir.join(sub).display()))?;
    }
    let record = RunRecord {
      id,
      title: title.unwrap_or_else(|| command.to_string()),
      command,
      config,
      status: RunStatus::Created,
      pinned: false,
      created_at: Utc::now(),
      started_at: None,
      finished_at: None,
      duration_seconds: None,
      treetime_version: version_info().version.to_owned(),
      inputs: vec![],
      config_hash: None,
      changed_settings: vec![],
      headline: BTreeMap::new(),
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
    let dir = self.run_dir(&record.id);
    let path = dir.join(RUN_FILE);
    let mut file =
      NamedTempFile::new_in(&dir).wrap_err_with(|| format!("When creating a temporary file in '{}'", dir.display()))?;
    file.write_all(format!("{}\n", json_write_str(record, JsonPretty(true))?).as_bytes())?;
    file
      .persist(&path)
      .map_err(|err| Report::new(err.error))
      .wrap_err_with(|| format!("When writing the run record '{}'", path.display()))?;
    Ok(())
  }

  pub fn list(&self) -> Result<Vec<RunRecord>, Report> {
    let records: Vec<RunRecord> = self.run_ids()?.iter().map(|id| self.read(id)).try_collect()?;
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

  pub fn trash(&self, id: &JobId) -> Result<(), Report> {
    self.read(id)?;
    let target = self.root.join(TRASH_DIR).join(id.as_str());
    if target.exists() {
      fs::remove_dir_all(&target).wrap_err_with(|| format!("When replacing the deleted run '{}'", target.display()))?;
    }
    fs::rename(self.run_dir(id), &target).wrap_err_with(|| format!("When deleting run `{}`", id.as_str()))
  }

  pub fn restore(&self, id: &JobId) -> Result<RunRecord, Report> {
    let source = self.root.join(TRASH_DIR).join(id.as_str());
    if !source.join(RUN_FILE).is_file() {
      return Err(not_found(format!("no deleted run with id `{}`", id.as_str())));
    }
    let target = self.run_dir(id);
    if target.exists() {
      return Err(conflict(format!("a run with id `{}` already exists", id.as_str())));
    }
    fs::rename(&source, &target).wrap_err_with(|| format!("When restoring run `{}`", id.as_str()))?;
    self.read(id)
  }

  pub fn purge(&self, id: &JobId) -> Result<(), Report> {
    let dir = self.root.join(TRASH_DIR).join(id.as_str());
    if !dir.is_dir() {
      return Err(not_found(format!("no deleted run with id `{}`", id.as_str())));
    }
    fs::remove_dir_all(&dir).wrap_err_with(|| format!("When purging run `{}`", id.as_str()))
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
      if name == TRASH_DIR || !entry.file_type()?.is_dir() || !entry.path().join(RUN_FILE).is_file() {
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
  let text = match fs::read_to_string(path) {
    Ok(text) => text,
    Err(err) if err.kind() == ErrorKind::NotFound => {
      return Err(not_found(format!("no run with id `{}`", id.as_str())));
    },
    Err(err) => return Err(Report::new(err).wrap_err(format!("When reading '{}'", path.display()))),
  };
  serde_json::from_str(&text).wrap_err_with(|| format!("When reading the run record '{}'", path.display()))
}
