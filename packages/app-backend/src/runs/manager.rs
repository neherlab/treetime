use crate::atomic_write::write_atomically;
use crate::command::{AppCommandExt, config_from_request, config_from_settings, ensure_in_apps};
use crate::config::catalog::command_settings;
use crate::job::{CancelToken, JobEvent, JobId, TerminalEvent};
use crate::runs::app_events::{AppChange, AppEventLog, run_stale_paths};
use crate::runs::downloads::{FileDownload, download_name};
use crate::runs::errors::{conflict, invalid};
use crate::runs::events::{EventLog, Subscriber, read_events};
use crate::runs::files::{RunFile, files_size, list_run_files, resolve_run_file};
use crate::runs::input_store::{
  InputStore, StoredInput, StoredUpload, UploadWriter, format_size, upload_too_large, validate_upload_name,
};
use crate::runs::lifecycle::StartedRun;
use crate::runs::package::{RunPackage, export_run};
use crate::runs::package_import::PackageImport;
use crate::runs::record::{
  CreateRunRequest, InputSet, RunInput, RunList, RunRecord, RunStatus, RunSummary, SaveRunRequest, StartRunRequest,
  UpdateRunRequest,
};
use crate::runs::run_name::{dataset_label, run_name};
use crate::runs::run_paths::{PortablePaths, resolve_settings};
use crate::runs::store::{NewRun, RunStore};
use app_commands::command_config::CommandConfig;
use app_commands::json_value::SparseConfig;
use chrono::{TimeDelta, Utc};
use deser::Serialize;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use parking_lot::{Mutex, MutexGuard};
use schemars::JsonSchema;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use treetime_utils::o;

const APP_EVENT_CAPACITY: usize = 1024;

const UPLOAD_AFTER_START: &str = "inputs can be uploaded only before the run starts";

const RELEASE_AFTER_START: &str = "its uploads can be removed only before the run starts";

const UNSTARTED_RUN_LIFETIME: TimeDelta = TimeDelta::days(1);

pub struct RunManager {
  store: RunStore,
  input_store: InputStore,
  active: Mutex<BTreeMap<JobId, Arc<ActiveRun>>>,
  pub(super) inputs: Mutex<BTreeMap<JobId, Arc<RunInputs>>>,
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
    store
      .remove_staged()
      .wrap_err("When removing the run packages that a previous process left unfinished")?;
    let input_store = store.input_store();
    input_store
      .remove_partial_uploads()
      .wrap_err("When removing the uploads that a previous process left unfinished")?;
    let first_app_seq =
      u64::try_from(Utc::now().timestamp_micros()).wrap_err("When numbering the app events from the current time")?;
    let manager = Arc::new(Self {
      store,
      input_store,
      active: Mutex::new(BTreeMap::new()),
      inputs: Mutex::new(BTreeMap::new()),
      records: Mutex::new(()),
      app_events: AppEventLog::new(APP_EVENT_CAPACITY, first_app_seq),
    });
    manager
      .remove_unreferenced_inputs()
      .wrap_err("When removing the uploaded inputs that no run references")?;
    Ok(manager)
  }

  pub fn app_events(&self) -> &AppEventLog {
    &self.app_events
  }

  pub fn store(&self) -> &RunStore {
    &self.store
  }

  pub fn input_store(&self) -> &InputStore {
    &self.input_store
  }

  pub fn list(&self) -> Result<RunList, Report> {
    Ok(RunList {
      runs: self
        .sized(self.store.list()?)?
        .iter()
        .filter(|record| record.config.command().in_apps())
        .map(RunRecord::summary)
        .collect(),
      active_runs: self.active_count(),
    })
  }

  pub fn active_count(&self) -> usize {
    self.active.lock().values().filter(|run| run.is_started()).count()
  }

  pub fn get(&self, id: &JobId) -> Result<RunRecord, Report> {
    self.with_size(self.store.read(id)?)
  }

  pub fn settings(&self, id: &JobId) -> Result<SparseConfig, Report> {
    let record = self.store.read(id)?;
    let mut settings = resolve_settings(&record, &self.store.run_dir(id), &self.input_store)?;
    record.config.command().remove_output_paths(&mut settings)?;
    Ok(SparseConfig(settings))
  }

  pub fn create(&self, request: &CreateRunRequest) -> Result<RunRecord, Report> {
    let id = JobId::random();
    let config = config_from_request(request.command, &request.config)?;
    let record = {
      let _records = self.records.lock();
      let PortableConfig {
        config,
        inputs,
        input_sets,
      } = self.portable(&id, &config)?;
      let dataset = dataset_label(request.example.as_deref(), &inputs);
      let title = run_name(dataset.as_deref(), &command_settings(request.command)?, &[]);
      self.store.create(NewRun {
        id,
        title,
        config,
        inputs,
        input_sets,
        source_run_id: request.source_run_id.clone(),
        example: request.example.clone(),
      })?
    };
    let record = self.with_size(record)?;
    self.app_events.append(
      AppChange::RunCreated { run: record.summary() },
      run_stale_paths(&record.id),
    );
    Ok(record)
  }

  pub fn start(
    self: &Arc<Self>,
    id: &JobId,
    request: &StartRunRequest,
    hook: ConfigHook,
  ) -> Result<StartedRun, Report> {
    let mut active = self.active.lock();
    let Some(inputs) = self.run_inputs(id)? else {
      return Err(conflict(format!("run `{}` has already started", id.as_str())));
    };
    let Some(_inputs) = inputs.idle() else {
      return Err(upload_in_progress(id, &inputs, "start"));
    };
    let record = self.modify(id, |record| {
      if record.status != RunStatus::Created {
        return Err(conflict(format!("run `{}` has already started", id.as_str())));
      }
      let command = request.command.unwrap_or_else(|| record.config.command());
      ensure_in_apps(command)?;
      let config = match &request.config {
        Some(config) => config_from_request(command, config)?,
        None if command != record.config.command() => {
          let current = resolve_settings(record, &self.store.run_dir(id), &self.input_store)?;
          config_from_request(command, &SparseConfig(current))?
        },
        None => record.config.clone(),
      };
      if request.config.is_some() || command != record.config.command() {
        let PortableConfig {
          config,
          inputs,
          input_sets,
        } = self.portable(id, &config)?;
        record.config = config;
        record.inputs = inputs;
        record.input_sets = input_sets;
      }
      if let Some(source) = &request.source_run_id {
        record.source_run_id = Some(source.clone());
      }
      if let Some(example) = &request.example {
        record.example = Some(example.clone());
      }
      record.uploads.clear();
      record.status = RunStatus::Running;
      record.started_at = Some(Utc::now());
      Ok(())
    })?;
    self.inputs.lock().remove(id);
    self.remove_unreferenced_inputs()?;
    let run = if let Some(run) = active.get(id) {
      Arc::clone(run)
    } else {
      let run = Arc::new(ActiveRun::open(&self.store.events_path(id))?);
      active.insert(id.clone(), Arc::clone(&run));
      run
    };
    run.started.store(true, Ordering::SeqCst);
    Ok(StartedRun::new(Arc::clone(self), record, run, hook))
  }

  pub fn cancel(&self, id: &JobId) -> Result<bool, Report> {
    let mut active = self.active.lock();
    if let Some(run) = active.get(id)
      && run.is_started()
    {
      run.token.cancel();
      return Ok(true);
    }
    let Some(inputs) = self.run_inputs(id)? else {
      return Ok(false);
    };
    let Some(_inputs) = inputs.idle() else {
      return Err(upload_in_progress(id, &inputs, "cancel"));
    };
    let run = match active.remove(id) {
      Some(run) => run,
      None => Arc::new(ActiveRun::open(&self.store.events_path(id))?),
    };
    run.events.append(JobEvent::Terminal {
      data: TerminalEvent::Cancelled { job_id: id.clone() },
    })?;
    self.modify(id, |record| {
      record.status = RunStatus::Cancelled;
      record.finished_at = Some(Utc::now());
      record.uploads.clear();
      Ok(())
    })?;
    self.inputs.lock().remove(id);
    self.remove_unreferenced_inputs()?;
    Ok(true)
  }

  pub fn delete(&self, id: &JobId) -> Result<(), Report> {
    let mut active = self.active.lock();
    let inputs = self.inputs.lock().get(id).cloned();
    let _upload = match &inputs {
      Some(inputs) => Some(inputs.idle().ok_or_else(|| upload_in_progress(id, inputs, "delete"))?),
      None => None,
    };
    {
      let _records = self.records.lock();
      let record = self.store.read(id)?;
      if record.status == RunStatus::Running || active.get(id).is_some_and(|run| run.is_started()) {
        return Err(conflict(format!(
          "run `{}` is computing; cancel it before deleting it",
          id.as_str()
        )));
      }
      self.store.remove(id)?;
    }
    active.remove(id);
    self.inputs.lock().remove(id);
    self.remove_unreferenced_inputs()?;
    self
      .app_events
      .append(AppChange::RunRemoved { id: id.clone() }, run_stale_paths(id));
    Ok(())
  }

  pub fn subscribe(&self, id: &JobId, from: u64, mut subscriber: Subscriber) -> Result<(), Report> {
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
        record.renamed = true;
      }
      if let Some(pinned) = request.pinned {
        record.pinned = pinned;
      }
      Ok(())
    })?;
    Ok(record.summary())
  }

  pub fn input_dirs(&self) -> Vec<PathBuf> {
    let root = self.input_store.root();
    if root.is_dir() {
      vec![root.to_path_buf()]
    } else {
      vec![]
    }
  }

  pub fn begin_upload(
    self: &Arc<Self>,
    id: &JobId,
    name: &str,
    replaces: &[PathBuf],
    max_total_size: usize,
  ) -> Result<InputUpload, Report> {
    validate_upload_name(name)?;
    let Some(inputs) = self.run_inputs(id)? else {
      return Err(started(id, UPLOAD_AFTER_START));
    };
    let _inputs = inputs.lock.lock();
    let record = self.store.read(id)?;
    if record.status != RunStatus::Created {
      return Err(started(id, UPLOAD_AFTER_START));
    }
    let replaced = self.replaced_uploads(id, &record, replaces)?;
    let retained = self.retained_size(&record, &replaced)?;
    let writer = self
      .input_store
      .begin(name, max_total_size.saturating_sub(retained), max_total_size)?;
    Ok(InputUpload {
      runs: Arc::clone(self),
      id: id.clone(),
      replaced,
      max_total_size,
      writer,
      uploading: Uploading::new(Arc::clone(&inputs), name),
    })
  }

  pub fn upload_allowance(
    &self,
    id: &JobId,
    replaces: &[PathBuf],
    max_total_size: usize,
  ) -> Result<UploadAllowance, Report> {
    let record = self.store.read(id)?;
    if record.status != RunStatus::Created {
      return Err(started(id, UPLOAD_AFTER_START));
    }
    let replaced = self.replaced_uploads(id, &record, replaces)?;
    let retained = self.retained_size(&record, &replaced)?;
    Ok(UploadAllowance {
      remaining: max_total_size.saturating_sub(retained),
      limit: max_total_size,
      limit_text: format_size(max_total_size),
    })
  }

  pub fn release_input(&self, id: &JobId, path: &Path) -> Result<(), Report> {
    let Some(inputs) = self.run_inputs(id)? else {
      return Err(started(id, RELEASE_AFTER_START));
    };
    let _inputs = inputs.lock.lock();
    let key = self.input_store.key_of(path).ok_or_else(|| not_an_upload(id, path))?;
    {
      let _records = self.records.lock();
      let mut record = self.store.read(id)?;
      if record.status != RunStatus::Created {
        return Err(started(id, RELEASE_AFTER_START));
      }
      record.uploads.retain(|upload| *upload != key);
      self.store.write(&record)?;
    }
    self.remove_unreferenced_inputs()
  }

  fn run_inputs(&self, id: &JobId) -> Result<Option<Arc<RunInputs>>, Report> {
    let mut inputs = self.inputs.lock();
    if let Some(run) = inputs.get(id) {
      return Ok(Some(Arc::clone(run)));
    }
    if self.store.read(id)?.status != RunStatus::Created {
      return Ok(None);
    }
    Ok(Some(Arc::clone(inputs.entry(id.clone()).or_default())))
  }

  fn replaced_uploads(
    &self,
    id: &JobId,
    record: &RunRecord,
    replaces: &[PathBuf],
  ) -> Result<BTreeSet<StoredInput>, Report> {
    replaces
      .iter()
      .map(|path| {
        self
          .input_store
          .key_of(path)
          .filter(|key| record.uploads.contains(key))
          .ok_or_else(|| not_an_upload(id, path))
      })
      .collect()
  }

  fn retained_size(&self, record: &RunRecord, replaced: &BTreeSet<StoredInput>) -> Result<usize, Report> {
    record
      .uploads
      .iter()
      .unique()
      .filter(|key| !replaced.contains(key))
      .map(|key| stored_size(&self.input_store, key))
      .sum()
  }

  fn commit_upload(&self, upload: StoredUpload, finished: &FinishedUpload) -> Result<StoredInput, Report> {
    let _records = self.records.lock();
    let mut record = self.store.read(finished.id)?;
    if record.status != RunStatus::Created {
      return Err(started(finished.id, UPLOAD_AFTER_START));
    }
    let retained = self.retained_size(&record, finished.replaced)?;
    if retained.saturating_add(upload.size()) > finished.max_total_size {
      return Err(upload_too_large(finished.max_total_size, finished.name));
    }
    let key = self.input_store.commit(upload)?;
    if !record.uploads.contains(&key) {
      record.uploads.push(key.clone());
      self.store.write(&record)?;
    }
    Ok(key)
  }

  pub fn files(&self, id: &JobId) -> Result<Vec<RunFile>, Report> {
    let record = self.store.read(id)?;
    list_run_files(&record, &self.store.out_dir(id))
  }

  pub fn file_path(&self, id: &JobId, relative: &str) -> Result<PathBuf, Report> {
    self.store.read(id)?;
    resolve_run_file(&self.store.out_dir(id), relative)
  }

  pub fn file_download(&self, id: &JobId, relative: &str) -> Result<FileDownload, Report> {
    let record = self.store.read(id)?;
    Ok(FileDownload {
      path: resolve_run_file(&self.store.out_dir(id), relative)?,
      name: download_name(&record, relative),
    })
  }

  pub fn save_file(&self, id: &JobId, relative: &str, destination: &Path) -> Result<(), Report> {
    let source = self.file_path(id, relative)?;
    write_atomically(destination, |file| {
      let mut reader = fs::File::open(&source).wrap_err_with(|| format!("When opening '{}'", source.display()))?;
      io::copy(&mut reader, file).wrap_err_with(|| format!("When copying '{}'", source.display()))?;
      Ok(())
    })
  }

  pub fn package(&self, id: &JobId) -> Result<RunPackage, Report> {
    export_run(self, id)
  }

  pub fn save_package(&self, id: &JobId, destination: &Path) -> Result<(), Report> {
    let package = self.package(id)?;
    write_atomically(destination, |file| package.write(file))
  }

  pub fn begin_import(self: &Arc<Self>, limit: usize) -> Result<PackageImport, Report> {
    PackageImport::begin(Arc::clone(self), limit)
  }

  pub fn publish_import(
    &self,
    mut record: RunRecord,
    staged: &Path,
    uploads: Vec<StoredUpload>,
    sets: &[(String, PathBuf)],
  ) -> Result<RunRecord, Report> {
    let published = {
      let _records = self.records.lock();
      if self.store.contains(&record.id) {
        let id = JobId::random();
        if record.group_id == record.id {
          record.group_id = id.clone();
        }
        record.id = id;
      }
      uploads
        .into_iter()
        .try_for_each(|upload| self.input_store.commit(upload).map(|_| ()))
        .and_then(|()| {
          sets
            .iter()
            .try_for_each(|(sha256, staged)| self.input_store.commit_set(staged, sha256))
        })
        .and_then(|()| self.store.publish(staged, &record))
    };
    if let Err(report) = published {
      self.remove_unreferenced_inputs()?;
      return Err(report);
    }
    let record = self.with_size(record)?;
    self.app_events.append(
      AppChange::RunCreated { run: record.summary() },
      run_stale_paths(&record.id),
    );
    Ok(record)
  }

  pub fn save(&self, id: &JobId, request: &SaveRunRequest) -> Result<(), Report> {
    if !request.destination.is_absolute() {
      return Err(invalid(format!(
        "the destination must be an absolute path, got '{}'",
        request.destination.display()
      )));
    }
    match &request.path {
      Some(relative) => self.save_file(id, relative, &request.destination),
      None => self.save_package(id, &request.destination),
    }
  }

  pub(super) fn modify(
    &self,
    id: &JobId,
    change: impl FnOnce(&mut RunRecord) -> Result<(), Report>,
  ) -> Result<RunRecord, Report> {
    let record = {
      let _records = self.records.lock();
      let mut record = self.store.read(id)?;
      change(&mut record)?;
      self.store.write(&record)?;
      record
    };
    let record = self.with_size(record)?;
    self
      .app_events
      .append(AppChange::RunUpdated { run: record.summary() }, run_stale_paths(id));
    Ok(record)
  }

  pub(super) fn finish(&self, id: &JobId, end: impl FnOnce()) {
    let mut active = self.active.lock();
    end();
    active.remove(id);
  }

  fn portable(&self, id: &JobId, config: &CommandConfig) -> Result<PortableConfig, Report> {
    let command = config.command();
    let mut settings = config.settings()?;
    let paths = PortablePaths::of_settings(command, &settings, &self.input_store)?.in_run(&self.store.run_dir(id));
    let inputs = paths.stored_inputs(command, &settings)?;
    paths.apply(command, &mut settings)?;
    Ok(PortableConfig {
      config: config_from_settings(command, &Value::Object(settings))?,
      inputs,
      input_sets: paths.input_sets(),
    })
  }

  fn with_size(&self, mut record: RunRecord) -> Result<RunRecord, Report> {
    let others = self.store.list()?;
    let references = others
      .iter()
      .filter(|other| other.id != record.id)
      .chain([&record])
      .flat_map(stored_references)
      .fold(BTreeMap::new(), count_reference);
    record.size = self.size_of(&record, &references)?;
    Ok(record)
  }

  fn sized(&self, records: Vec<RunRecord>) -> Result<Vec<RunRecord>, Report> {
    let references = records
      .iter()
      .flat_map(stored_references)
      .fold(BTreeMap::new(), count_reference);
    records
      .into_iter()
      .map(|mut record| {
        record.size = self.size_of(&record, &references)?;
        Ok(record)
      })
      .collect()
  }

  fn size_of(&self, record: &RunRecord, references: &BTreeMap<StoredEntry, usize>) -> Result<usize, Report> {
    let outputs = files_size(&self.store.out_dir(&record.id))?;
    let own_inputs = stored_references(record)
      .iter()
      .filter(|entry| references.get(*entry) == Some(&1))
      .map(|entry| self.entry_size(entry))
      .sum::<Result<usize, _>>()?;
    Ok(outputs + own_inputs)
  }

  fn entry_size(&self, entry: &StoredEntry) -> Result<usize, Report> {
    match entry {
      StoredEntry::File(key) => file_size(&self.input_store.path(key)),
      StoredEntry::Set(sha256) => self.input_store.set_size(sha256),
    }
  }

  fn remove_unreferenced_inputs(&self) -> Result<(), Report> {
    let _records = self.records.lock();
    let referenced = self
      .store
      .list()?
      .iter()
      .flat_map(stored_references)
      .collect::<BTreeSet<_>>();
    for key in self.input_store.keys()? {
      if !referenced.contains(&StoredEntry::File(key.clone())) {
        self.input_store.remove(&key)?;
      }
    }
    for sha256 in self.input_store.set_keys()? {
      if !referenced.contains(&StoredEntry::Set(sha256.clone())) {
        self.input_store.remove_set(&sha256)?;
      }
    }
    Ok(())
  }
}

pub type ConfigHook = Box<dyn Fn(&mut Value) -> Result<(), Report> + Send>;

/// A file uploaded for a run that has not started.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct UploadedInput {
  /// File name of the upload.
  pub name: String,
  /// Path to use for the file in the configurations of runs, checks, and other requests.
  pub path: PathBuf,
  /// Size in bytes.
  pub size: usize,
  /// SHA-256 of the contents, as lowercase hexadecimal.
  pub sha256: String,
}

/// Bytes that the uploads of a run that has not started may add before its inputs reach the upload limit.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct UploadAllowance {
  /// Bytes that the next uploads may add in total; the files that the request names as replaced do not count.
  pub remaining: usize,
  /// Limit of the total size of the inputs of a run, in bytes.
  pub limit: usize,
  /// The limit as the upload errors name it, such as `1 GiB`.
  pub limit_text: String,
}

pub struct InputUpload {
  runs: Arc<RunManager>,
  id: JobId,
  replaced: BTreeSet<StoredInput>,
  max_total_size: usize,
  writer: UploadWriter,
  uploading: Uploading,
}

impl InputUpload {
  pub fn write(&mut self, chunk: &[u8]) -> Result<(), Report> {
    self.writer.write(chunk)
  }

  pub fn finish(self) -> Result<UploadedInput, Report> {
    let Self {
      runs,
      id,
      replaced,
      max_total_size,
      writer,
      uploading,
    } = self;
    let upload = writer.finish()?;
    let size = upload.size();
    let key = {
      let _inputs = uploading.inputs.lock.lock();
      runs.commit_upload(
        upload,
        &FinishedUpload {
          id: &id,
          name: &uploading.name,
          replaced: &replaced,
          max_total_size,
        },
      )?
    };
    Ok(UploadedInput {
      name: uploading.name.clone(),
      path: runs.input_store.path(&key),
      size,
      sha256: key.sha256,
    })
  }
}

pub(super) struct ActiveRun {
  pub(super) token: CancelToken,
  pub(super) events: EventLog,
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

#[derive(Default)]
pub(super) struct RunInputs {
  lock: Mutex<()>,
  uploading: Mutex<Vec<String>>,
}

impl RunInputs {
  fn idle(&self) -> Option<MutexGuard<'_, ()>> {
    let guard = self.lock.try_lock()?;
    self.uploading.lock().is_empty().then_some(guard)
  }
}

struct Uploading {
  inputs: Arc<RunInputs>,
  name: String,
}

impl Uploading {
  fn new(inputs: Arc<RunInputs>, name: &str) -> Self {
    inputs.uploading.lock().push(name.to_owned());
    Self {
      inputs,
      name: name.to_owned(),
    }
  }
}

impl Drop for Uploading {
  fn drop(&mut self) {
    let mut uploading = self.inputs.uploading.lock();
    if let Some(index) = uploading.iter().position(|name| *name == self.name) {
      uploading.remove(index);
    }
  }
}

struct FinishedUpload<'a> {
  id: &'a JobId,
  name: &'a str,
  replaced: &'a BTreeSet<StoredInput>,
  max_total_size: usize,
}

struct PortableConfig {
  config: CommandConfig,
  inputs: Vec<RunInput>,
  input_sets: Vec<InputSet>,
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum StoredEntry {
  File(StoredInput),
  Set(String),
}

fn stored_references(record: &RunRecord) -> BTreeSet<StoredEntry> {
  record
    .uploads
    .iter()
    .cloned()
    .chain(record.inputs.iter().filter_map(RunInput::stored))
    .map(StoredEntry::File)
    .chain(record.input_sets.iter().map(|set| StoredEntry::Set(set.sha256.clone())))
    .collect()
}

fn count_reference(mut counts: BTreeMap<StoredEntry, usize>, entry: StoredEntry) -> BTreeMap<StoredEntry, usize> {
  *counts.entry(entry).or_default() += 1;
  counts
}

fn stored_size(store: &InputStore, key: &StoredInput) -> Result<usize, Report> {
  file_size(&store.path(key))
}

fn upload_in_progress(id: &JobId, inputs: &RunInputs, action: &str) -> Report {
  let uploading = inputs.uploading.lock().first().cloned();
  conflict(format!(
    "run `{}` cannot {action} while {} is uploaded; {action} it again when the upload ends",
    id.as_str(),
    uploading.map_or_else(|| o!("a file"), |name| format!("`{name}`"))
  ))
}

fn started(id: &JobId, consequence: &str) -> Report {
  conflict(format!("run `{}` has already started; {consequence}", id.as_str()))
}

fn not_an_upload(id: &JobId, path: &Path) -> Report {
  invalid(format!(
    "'{}' is not an upload of run `{}`",
    path.display(),
    id.as_str()
  ))
}

fn file_size(path: &Path) -> Result<usize, Report> {
  match fs::metadata(path) {
    Ok(metadata) if metadata.is_file() => Ok(usize::try_from(metadata.len())?),
    Ok(_) => Ok(0),
    Err(err) if err.kind() == io::ErrorKind::NotFound => Ok(0),
    Err(err) => Err(Report::new(err).wrap_err(format!("When reading the size of '{}'", path.display()))),
  }
}
