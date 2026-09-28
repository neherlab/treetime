use crate::command::{CommandOutcome, OutputFile};
use crate::job::{CancelToken, JobEvent, JobId, JobProgress, JobStarted, TerminalEvent, run_job};
use crate::runs::app_events::{AppChange, AppEventLog, run_stale_paths};
use crate::runs::errors::{UploadTooLarge, conflict, invalid};
use crate::runs::events::{EventLog, Subscriber, read_events};
use crate::runs::files::{RunFile, list_run_files, resolve_run_file, write_run_zip};
use crate::runs::headline::{RunHeadline, run_headline};
use crate::runs::inputs::{file_sha256, hash_inputs};
use crate::runs::record::{
  CreateRunRequest, RunError, RunList, RunRecord, RunStatus, RunSummary, StartRunRequest, UpdateRunRequest,
};
use crate::runs::store::RunStore;
use chrono::{TimeDelta, Utc};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use log::error;
use parking_lot::Mutex;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::BTreeMap;
use std::fs;
use std::io::{Cursor, Read, Write};
use std::path::{Component, Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;
use tempfile::NamedTempFile;
use treetime::cancel::Cancel;
use treetime_utils::fmt::float::float_to_significant_digits;

const UPLOAD_BUFFER_SIZE: usize = 1 << 16;

const APP_EVENT_CAPACITY: usize = 1024;

const UNSTARTED_RUN_LIFETIME: TimeDelta = TimeDelta::days(1);

pub type ConfigHook = Box<dyn FnOnce(&mut Value) -> Result<(), Report> + Send>;

pub fn unconfined() -> ConfigHook {
  Box::new(|_config: &mut Value| Ok(()))
}

pub struct RunManager {
  store: RunStore,
  active: Mutex<BTreeMap<JobId, Arc<ActiveRun>>>,
  records: Mutex<()>,
  app_events: AppEventLog,
}

impl RunManager {
  pub fn open(root: &Path) -> Result<Arc<Self>, Report> {
    let store = RunStore::open(root)?;
    store
      .recover_interrupted()
      .wrap_err("When marking the runs of a previous process as interrupted")?;
    store
      .remove_unstarted(Utc::now() - UNSTARTED_RUN_LIFETIME)
      .wrap_err("When removing runs that were created for uploads and never started")?;
    let first_app_seq =
      usize::try_from(Utc::now().timestamp_micros()).wrap_err("When numbering the app events from the current time")?;
    Ok(Arc::new(Self {
      store,
      active: Mutex::new(BTreeMap::new()),
      records: Mutex::new(()),
      app_events: AppEventLog::new(APP_EVENT_CAPACITY, first_app_seq),
    }))
  }

  pub fn app_events(&self) -> &AppEventLog {
    &self.app_events
  }

  pub fn store(&self) -> &RunStore {
    &self.store
  }

  pub fn list(&self) -> Result<RunList, Report> {
    Ok(RunList {
      runs: self.store.list()?.iter().map(RunRecord::summary).collect(),
      active_runs: self.active_count(),
    })
  }

  pub fn active_count(&self) -> usize {
    self.active.lock().values().filter(|run| run.is_started()).count()
  }

  pub fn get(&self, id: &JobId) -> Result<RunRecord, Report> {
    self.store.read(id)
  }

  pub fn create(&self, request: CreateRunRequest) -> Result<RunRecord, Report> {
    let record = self.store.create(request.command, request.config)?;
    self.app_events.append(
      AppChange::RunCreated { run: record.summary() },
      run_stale_paths(&record.id),
    );
    Ok(record)
  }

  pub fn start(self: &Arc<Self>, id: &JobId, request: StartRunRequest, hook: ConfigHook) -> Result<StartedRun, Report> {
    let mut active = self.active.lock();
    let record = self.modify(id, |record| {
      if record.status != RunStatus::Created {
        return Err(conflict(format!("run `{}` has already started", id.as_str())));
      }
      if let Some(command) = request.command {
        record.command = command;
      }
      if let Some(config) = request.config {
        let Value::Object(config) = config else {
          return Err(invalid("a command configuration must be a mapping of settings"));
        };
        record.config = config;
      }
      record.status = RunStatus::Running;
      record.started_at = Some(Utc::now());
      Ok(())
    })?;
    let run = if let Some(run) = active.get(id) {
      Arc::clone(run)
    } else {
      let run = Arc::new(ActiveRun::open(&self.store.events_path(id))?);
      active.insert(id.clone(), Arc::clone(&run));
      run
    };
    run.started.store(true, Ordering::SeqCst);
    Ok(StartedRun {
      manager: Arc::clone(self),
      record,
      run,
      hook,
    })
  }

  pub fn cancel(&self, id: &JobId) -> Result<bool, Report> {
    let mut active = self.active.lock();
    if let Some(run) = active.get(id)
      && run.is_started()
    {
      run.token.cancel();
      return Ok(true);
    }
    let record = self.store.read(id)?;
    if record.status != RunStatus::Created {
      return Ok(false);
    }
    let run = match active.remove(id) {
      Some(run) => run,
      None => Arc::new(ActiveRun::open(&self.store.events_path(id))?),
    };
    run
      .events
      .append(JobEvent::Terminal(TerminalEvent::Cancelled { job_id: id.clone() }))?;
    self.modify(id, |record| {
      record.status = RunStatus::Cancelled;
      record.finished_at = Some(Utc::now());
      Ok(())
    })?;
    Ok(true)
  }

  pub fn subscribe(&self, id: &JobId, from: usize, mut subscriber: Subscriber) -> Result<(), Report> {
    let mut active = self.active.lock();
    if let Some(run) = active.get(id) {
      run.events.subscribe(from, subscriber);
      return Ok(());
    }
    let record = self.store.read(id)?;
    if record.status == RunStatus::Created {
      let run = Arc::new(ActiveRun::open(&self.store.events_path(id))?);
      run.events.subscribe(from, subscriber);
      active.insert(id.clone(), run);
      return Ok(());
    }
    for event in read_events(&self.store.events_path(id), from)? {
      if !subscriber(&event) {
        break;
      }
    }
    Ok(())
  }

  pub fn update(&self, id: &JobId, request: UpdateRunRequest) -> Result<RunSummary, Report> {
    let record = self.modify(id, |record| {
      if let Some(title) = request.title {
        if title.trim().is_empty() {
          return Err(invalid("a run title must not be empty"));
        }
        record.title = title;
      }
      if let Some(pinned) = request.pinned {
        record.pinned = pinned;
      }
      Ok(())
    })?;
    Ok(record.summary())
  }

  pub fn input_dirs(&self) -> Result<Vec<PathBuf>, Report> {
    self.store.input_dirs()
  }

  pub fn upload_input(
    &self,
    id: &JobId,
    name: &str,
    body: &mut dyn Read,
    max_total_size: usize,
  ) -> Result<UploadedInput, Report> {
    let record = self.store.read(id)?;
    if record.status != RunStatus::Created {
      return Err(conflict(format!(
        "run `{}` has already started; inputs can be uploaded only before the run starts",
        id.as_str()
      )));
    }
    validate_upload_name(name)?;
    let dir = self.store.inputs_dir(id);
    let target = dir.join(name);
    let others: usize = fs::read_dir(&dir)
      .wrap_err_with(|| format!("When listing '{}'", dir.display()))?
      .map(|entry| -> Result<usize, Report> {
        let entry = entry?;
        let metadata = entry.metadata()?;
        Ok(if metadata.is_file() && entry.path() != target {
          usize::try_from(metadata.len())?
        } else {
          0
        })
      })
      .sum::<Result<usize, Report>>()?;
    let allowed = max_total_size.saturating_sub(others);

    let mut file =
      NamedTempFile::new_in(&dir).wrap_err_with(|| format!("When creating a file in '{}'", dir.display()))?;
    let mut buffer = vec![0_u8; UPLOAD_BUFFER_SIZE];
    let mut size: usize = 0;
    loop {
      let read = body.read(&mut buffer).wrap_err("When receiving the upload")?;
      if read == 0 {
        break;
      }
      size += read;
      if size > allowed {
        return Err(Report::new(UploadTooLarge {
          limit: format_size(max_total_size),
          name: name.to_owned(),
        }));
      }
      file.write_all(&buffer[..read]).wrap_err("When storing the upload")?;
    }
    file
      .persist(&target)
      .map_err(|err| Report::new(err.error))
      .wrap_err_with(|| format!("When storing the upload as '{}'", target.display()))?;
    let (size, sha256) = file_sha256(&target)?;
    Ok(UploadedInput {
      name: name.to_owned(),
      path: target,
      size,
      sha256,
    })
  }

  pub fn files(&self, id: &JobId) -> Result<Vec<RunFile>, Report> {
    let record = self.store.read(id)?;
    list_run_files(&self.store.out_dir(id), &record.output_files)
  }

  pub fn file_path(&self, id: &JobId, relative: &str) -> Result<PathBuf, Report> {
    self.store.read(id)?;
    resolve_run_file(&self.store.out_dir(id), relative)
  }

  pub fn zip(&self, id: &JobId) -> Result<Vec<u8>, Report> {
    self.store.read(id)?;
    let mut buffer = Cursor::new(vec![]);
    write_run_zip(&self.store.out_dir(id), id.as_str(), &mut buffer)?;
    Ok(buffer.into_inner())
  }

  fn modify(&self, id: &JobId, change: impl FnOnce(&mut RunRecord) -> Result<(), Report>) -> Result<RunRecord, Report> {
    let _records = self.records.lock();
    let mut record = self.store.read(id)?;
    change(&mut record)?;
    self.store.write(&record)?;
    self
      .app_events
      .append(AppChange::RunUpdated { run: record.summary() }, run_stale_paths(id));
    Ok(record)
  }

  fn finish(&self, id: &JobId) {
    self.active.lock().remove(id);
  }
}

/// A file uploaded into a run's `inputs/` folder.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct UploadedInput {
  /// File name inside the run's `inputs/` folder.
  pub name: String,
  /// Path to use for the file in the run's configuration.
  pub path: PathBuf,
  /// Size in bytes.
  pub size: usize,
  /// SHA-256 of the contents, as lowercase hexadecimal.
  pub sha256: String,
}

pub struct StartedRun {
  manager: Arc<RunManager>,
  record: RunRecord,
  run: Arc<ActiveRun>,
  hook: ConfigHook,
}

impl StartedRun {
  pub fn record(&self) -> &RunRecord {
    &self.record
  }

  pub fn run(self) -> TerminalEvent {
    let Self {
      manager,
      record,
      run,
      hook,
    } = self;
    let id = record.id.clone();
    let command = record.command;
    let out_dir = manager.store.out_dir(&id);
    let emit = |event: JobEvent| {
      if let Err(err) = run.events.append(event) {
        error!("When recording an event of run {}: {err:#}", id.as_str());
      }
    };
    emit(JobEvent::Started(JobStarted {
      job_id: id.clone(),
      command,
    }));
    let progress = JobProgress::new(emit);
    let clock = Instant::now();

    let terminal = run_job(&id, &run.token, || {
      let mut config = Value::Object(record.config.clone());
      hook(&mut config)?;
      let prepared = command.prepare_run(&config, &out_dir)?;
      let hashed = hash_inputs(command, &prepared.config)?;
      manager.modify(&id, |record| {
        record.config = prepared.config.clone();
        record.changed_settings = prepared.changed_settings.clone();
        record.inputs = hashed.inputs;
        record.config_hash = Some(hashed.config_hash);
        Ok(())
      })?;
      run.token.check()?;
      prepared.args.run(&run.token, &progress)
    });

    let duration = clock.elapsed().as_secs_f64();
    let finished = manager.modify(&id, |record| {
      record.finished_at = Some(Utc::now());
      record.duration_seconds = Some(duration);
      match &terminal {
        TerminalEvent::Ok { result, .. } => {
          record.status = RunStatus::Ok;
          record.output_files = relative_outputs(result, &out_dir);
          record.headline = run_headline(record, &out_dir).unwrap_or_else(|err| {
            error!(
              "When reading the headline results of run {}: {err:#}",
              record.id.as_str()
            );
            RunHeadline::default()
          });
        },
        TerminalEvent::Error { message, causes, .. } => {
          record.status = RunStatus::Error;
          record.error = Some(RunError {
            message: message.clone(),
            causes: causes.clone(),
          });
        },
        TerminalEvent::Cancelled { .. } => record.status = RunStatus::Cancelled,
        TerminalEvent::Interrupted { .. } => record.status = RunStatus::Interrupted,
      }
      Ok(())
    });
    if let Err(err) = finished {
      error!("When recording the end of run {}: {err:#}", id.as_str());
    }
    emit(JobEvent::Terminal(terminal.clone()));
    manager.finish(&id);
    terminal
  }
}

struct ActiveRun {
  token: CancelToken,
  events: EventLog,
  started: AtomicBool,
}

impl ActiveRun {
  fn open(events_path: &Path) -> Result<Self, Report> {
    Ok(Self {
      token: CancelToken::default(),
      events: EventLog::open(events_path)?,
      started: AtomicBool::new(false),
    })
  }

  fn is_started(&self) -> bool {
    self.started.load(Ordering::SeqCst)
  }
}

fn relative_outputs(outcome: &CommandOutcome, out_dir: &Path) -> Vec<OutputFile> {
  outcome
    .output_files
    .iter()
    .map(|file| OutputFile {
      path: file.path.strip_prefix(out_dir).unwrap_or(&file.path).to_path_buf(),
      kind: file.kind,
    })
    .collect_vec()
}

fn validate_upload_name(name: &str) -> Result<(), Report> {
  let mut components = Path::new(name).components();
  let single = matches!(components.next(), Some(Component::Normal(_))) && components.next().is_none();
  if !single || name.starts_with('.') || name.contains(['/', '\\']) {
    return Err(invalid(format!(
      "upload name `{name}` must be a plain file name without directories"
    )));
  }
  Ok(())
}

#[allow(
  clippy::as_conversions,
  clippy::cast_precision_loss,
  reason = "sizes are formatted for a message; precision beyond two decimals is irrelevant"
)]
fn format_size(bytes: usize) -> String {
  const UNITS: [&str; 5] = ["bytes", "KiB", "MiB", "GiB", "TiB"];
  let (value, unit) = UNITS
    .iter()
    .skip(1)
    .fold((bytes as f64, UNITS[0]), |(value, unit), next| {
      if value >= 1024.0 {
        (value / 1024.0, *next)
      } else {
        (value, unit)
      }
    });
  if unit == UNITS[0] {
    format!("{bytes} bytes")
  } else {
    format!("{} {unit}", float_to_significant_digits(value, 3))
  }
}
