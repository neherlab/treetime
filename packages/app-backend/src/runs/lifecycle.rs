use crate::command::AppCommandExt;
use crate::command::config_from_settings;
use crate::command::run_command;
use crate::command::{CommandOutcome, OutputAlignment, OutputFile};
use crate::config::catalog::{OUTPUT_SELECTION_KEY, command_settings};
use crate::job::{JobEvent, JobId, JobProgress, JobStarted, TerminalEvent, run_job, run_job_step};
use crate::run_checks::named_settings;
use crate::runs::errors::invalid;
use crate::runs::headline::{RunHeadline, run_headline};
use crate::runs::inputs::hash_inputs;
use crate::runs::manager::{ActiveRun, ConfigHook, RunManager};
use crate::runs::record::{RunError, RunInput, RunRecord, RunStatus};
use crate::runs::run_name::{dataset_label, run_name};
use crate::runs::run_paths::{OUT_DIR, PortablePaths, resolve_settings};
use crate::runs::setting_differences::{ComparedConfig, compared_inputs, setting_differences};
use crate::runs::warnings::WarningCollector;
use app_commands::command::AppCommand;
use app_commands::json_value::SparseConfig;
use chrono::Utc;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use log::{error, info};
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::iter;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use treetime::cancel::Cancel;
use treetime::progress::RunWarning;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::make_error;

const NO_ERROR_LINES: &str = "the run worker stopped without an error message";

const SEED_KEY: &str = "seed";

const MAX_DRAWN_SEED: u64 = (1 << 53) - 1;

pub fn compute_run(
  id: &JobId,
  command: AppCommand,
  config: &Map<String, Value>,
  out_dir: &Path,
  cancel: &dyn Cancel,
  emit: impl Fn(JobEvent) + Send + Sync,
) -> RunOutcome {
  let progress = JobProgress::new(emit);
  let log = WarningCollector::new(&progress);
  let terminal = run_job(id, cancel, || {
    let prepared = command.prepare_run(&Value::Object(config.clone()), out_dir)?;
    run_command(&prepared.args, cancel, &progress, &log)
  });
  RunOutcome {
    terminal,
    warnings: log.into_warnings(),
  }
}

pub struct StartedRun {
  manager: Arc<RunManager>,
  record: RunRecord,
  run: Arc<ActiveRun>,
  hook: ConfigHook,
}

impl StartedRun {
  pub(super) fn new(manager: Arc<RunManager>, record: RunRecord, run: Arc<ActiveRun>, hook: ConfigHook) -> Self {
    Self {
      manager,
      record,
      run,
      hook,
    }
  }

  pub fn record(&self) -> &RunRecord {
    &self.record
  }

  pub fn run(self) -> TerminalEvent {
    match self.prepare() {
      Prepared::Ready(prepared) => {
        let outcome = compute_run(
          &prepared.job.id,
          prepared.job.command,
          &prepared.local_config,
          &prepared.job.out_dir(),
          &prepared.run.token,
          |event| prepared.emit(event),
        );
        prepared.complete(outcome)
      },
      Prepared::Ended(terminal) => terminal,
    }
  }

  pub fn prepare(self) -> Prepared {
    let Self {
      manager,
      record,
      run,
      hook,
    } = self;
    let clock = Instant::now();
    let id = record.id.clone();
    let command = record.config.command();
    let run_dir = manager.store().run_dir(&id);
    append_event(
      &run,
      &id,
      JobEvent::Started {
        data: JobStarted {
          job_id: id.clone(),
          command,
        },
      },
    );
    let job = run_job_step(&id, &run.token, || {
      let store = manager.input_store();
      let mut config = Value::Object(resolve_settings(
        command,
        &record.config.settings()?,
        &record.inputs,
        &run_dir,
        store,
      )?);
      hook(&mut config)?;
      let Value::Object(settings) = &config else {
        return make_error!("a command configuration must be a mapping of settings");
      };
      let paths = PortablePaths::of_settings(command, settings, store)?.in_run(&run_dir);
      let prepared = command.prepare_run(&config, &run_dir.join(OUT_DIR))?;
      let hashed = hash_inputs(command, &prepared.config)?;
      let requested = requested_outputs(&config);
      let seed = prepared.args.random_step().map(|_| {
        prepared
          .args
          .seed()
          .unwrap_or_else(|| rand::random_range(0..=MAX_DRAWN_SEED))
      });
      let mut local_config = prepared.config.clone();
      if let Some(seed) = seed {
        local_config.insert(SEED_KEY.to_owned(), Value::from(seed));
      }
      let mut portable = local_config.clone();
      paths.apply(command, &mut portable)?;
      let resolved = config_from_settings(command, &Value::Object(portable.clone()))?;
      let sources = hashed
        .inputs
        .iter()
        .map(|input| (paths.portable(&input.path), input.path.clone()))
        .collect::<BTreeMap<_, _>>();
      let inputs = hashed
        .inputs
        .into_iter()
        .map(|input| RunInput {
          path: paths.portable(&input.path),
          ..input
        })
        .collect_vec();
      let group_id = derivation_group(&manager, &record, &inputs)?;
      let title = run_title(&manager, &record, &prepared.config, &inputs, &hook)?;
      manager.modify(&id, |record| {
        if !record.renamed {
          record.title = title;
        }
        record.config = resolved;
        record.changed_settings = prepared
          .changed_settings
          .iter()
          .filter(|key| requested.is_some() || key.as_str() != OUTPUT_SELECTION_KEY)
          .cloned()
          .collect();
        record.inputs = inputs.clone();
        record.config_hash = Some(hashed.config_hash);
        record.seed = seed;
        record.group_id = group_id;
        Ok(())
      })?;
      run.token.check()?;
      Ok(PreparedJob {
        job: RunJob {
          id: id.clone(),
          command,
          config: SparseConfig(portable),
          run_dir: run_dir.clone(),
          inputs: job_inputs(&inputs),
        },
        local_config,
        sources,
        paths,
      })
    });
    match job {
      Ok(PreparedJob {
        job,
        local_config,
        sources,
        paths,
      }) => Prepared::Ready(Box::new(PreparedRun {
        manager,
        run,
        clock,
        job,
        local_config,
        sources,
        paths,
      })),
      Err(terminal) => {
        let ended = RunEnd {
          manager: &manager,
          run: &run,
          id: &id,
          run_dir: &run_dir,
          paths: None,
          clock,
        };
        Prepared::Ended(ended.record(RunOutcome {
          terminal,
          warnings: vec![],
        }))
      },
    }
  }
}

pub enum Prepared {
  Ready(Box<PreparedRun>),
  Ended(TerminalEvent),
}

pub struct PreparedRun {
  manager: Arc<RunManager>,
  run: Arc<ActiveRun>,
  clock: Instant,
  job: RunJob,
  local_config: Map<String, Value>,
  sources: BTreeMap<PathBuf, PathBuf>,
  paths: PortablePaths,
}

impl PreparedRun {
  pub fn job(&self) -> &RunJob {
    &self.job
  }

  pub fn input_source(&self, path: &Path) -> Option<&Path> {
    self.sources.get(path).map(PathBuf::as_path)
  }

  pub fn is_cancelled(&self) -> bool {
    self.run.token.is_cancelled()
  }

  pub fn append_job_event(&self, event: JobEvent) -> Result<(), Report> {
    if matches!(event, JobEvent::Started { .. } | JobEvent::Terminal { .. }) {
      let kind: &'static str = (&event).into();
      return Err(invalid(format!(
        "a run worker cannot report a `{kind}` event of run `{}`; only the back end starts and ends a run",
        self.job.id.as_str()
      )));
    }
    self.run.events.append(event)?;
    Ok(())
  }

  pub fn complete(self, outcome: RunOutcome) -> TerminalEvent {
    let ended = RunEnd {
      manager: &self.manager,
      run: &self.run,
      id: &self.job.id,
      run_dir: &self.job.run_dir,
      paths: Some(&self.paths),
      clock: self.clock,
    };
    ended.record(outcome)
  }

  pub fn fail(self, lines: Vec<String>) -> TerminalEvent {
    let job_id = self.job.id.clone();
    let terminal = if self.is_cancelled() {
      TerminalEvent::Cancelled { job_id }
    } else {
      let mut causes = lines;
      let message = causes.pop().unwrap_or_else(|| NO_ERROR_LINES.to_owned());
      TerminalEvent::Error {
        job_id,
        message,
        causes,
      }
    };
    self.complete(RunOutcome {
      terminal,
      warnings: vec![],
    })
  }

  fn emit(&self, event: JobEvent) {
    append_event(&self.run, &self.job.id, event);
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RunJob {
  pub id: JobId,
  pub command: AppCommand,
  pub config: SparseConfig,
  pub run_dir: PathBuf,
  pub inputs: Vec<JobInput>,
}

impl RunJob {
  pub fn out_dir(&self) -> PathBuf {
    self.run_dir.join(OUT_DIR)
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct JobInput {
  pub path: PathBuf,
  pub size: usize,
  pub sha256: String,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RunOutcome {
  pub terminal: TerminalEvent,
  pub warnings: Vec<RunWarning>,
}

struct PreparedJob {
  job: RunJob,
  local_config: Map<String, Value>,
  sources: BTreeMap<PathBuf, PathBuf>,
  paths: PortablePaths,
}

struct RunEnd<'a> {
  manager: &'a RunManager,
  run: &'a ActiveRun,
  id: &'a JobId,
  run_dir: &'a Path,
  paths: Option<&'a PortablePaths>,
  clock: Instant,
}

impl RunEnd<'_> {
  fn record(&self, outcome: RunOutcome) -> TerminalEvent {
    let RunOutcome { terminal, warnings } = outcome;
    let duration = self.clock.elapsed().as_secs_f64();
    self.manager.finish(self.id, || {
      let finished = self.manager.modify(self.id, |record| {
        record.finished_at = Some(Utc::now());
        record.warnings = warnings;
        record.duration_seconds = Some(duration);
        match &terminal {
          TerminalEvent::Ok { result, .. } => {
            record.status = RunStatus::Ok;
            record.output_files = relative_outputs(result, self.run_dir, self.paths);
            record.headline = run_headline(record, self.run_dir, self.manager.input_store()).unwrap_or_else(|err| {
              error!(
                "When reading the headline results of run {}: {err:#}",
                record.id.as_str()
              );
              RunHeadline::default()
            });
          },
          TerminalEvent::Error { message, causes, .. } => {
            record.status = RunStatus::Error;
            let catalog = command_settings(record.config.command())?;
            record.error = Some(RunError {
              message: message.clone(),
              causes: causes.clone(),
              settings: named_settings(iter::once(message).chain(causes).map(String::as_str), &catalog.settings),
            });
          },
          TerminalEvent::Cancelled { .. } => record.status = RunStatus::Cancelled,
          TerminalEvent::Interrupted { .. } => record.status = RunStatus::Interrupted,
        }
        Ok(())
      });
      if let Err(err) = finished {
        error!("When recording the end of run {}: {err:#}", self.id.as_str());
      }
      append_event(self.run, self.id, JobEvent::Terminal { data: terminal.clone() });
    });
    match json_write_str(&terminal, JsonPretty(false)) {
      Ok(text) => info!("Run {} ended: {text}", self.id.as_str()),
      Err(err) => error!(
        "Run {} ended; its terminal event cannot be shown: {err}",
        self.id.as_str()
      ),
    }
    terminal
  }
}

fn append_event(run: &ActiveRun, id: &JobId, event: JobEvent) {
  if let Err(err) = run.events.append(event) {
    error!("When recording an event of run {}: {err:#}", id.as_str());
  }
}

fn job_inputs(inputs: &[RunInput]) -> Vec<JobInput> {
  inputs
    .iter()
    .unique_by(|input| &input.path)
    .map(|input| JobInput {
      path: input.path.clone(),
      size: input.size,
      sha256: input.sha256.clone(),
    })
    .collect_vec()
}

fn relative_outputs(outcome: &CommandOutcome, run_dir: &Path, paths: Option<&PortablePaths>) -> Vec<OutputFile> {
  outcome
    .output_files
    .iter()
    .map(|file| {
      let path = file.path.strip_prefix(run_dir).unwrap_or(&file.path);
      OutputFile {
        path: path.strip_prefix(OUT_DIR).unwrap_or(path).to_path_buf(),
        kind: file.kind,
        alignment: file.alignment.as_ref().map(|alignment| OutputAlignment {
          number: alignment.number,
          path: paths.map_or_else(|| alignment.path.clone(), |paths| paths.portable(&alignment.path)),
        }),
      }
    })
    .collect_vec()
}

fn derivation_group(manager: &RunManager, record: &RunRecord, inputs: &[RunInput]) -> Result<JobId, Report> {
  let records = manager.store().list()?;
  if let Some(source) = record
    .source_run_id
    .as_ref()
    .and_then(|source| records.iter().find(|other| &other.id == source))
  {
    return Ok(source.group_id.clone());
  }
  let digests = input_digests(inputs);
  if digests.is_empty() {
    return Ok(record.id.clone());
  }
  Ok(
    records
      .iter()
      .filter(|other| other.id != record.id && other.created_at < record.created_at && other.config_hash.is_some())
      .filter(|other| input_digests(&other.inputs) == digests)
      .min_by_key(|other| (other.created_at, other.id.clone()))
      .map_or_else(|| record.id.clone(), |other| other.group_id.clone()),
  )
}

fn run_title(
  manager: &RunManager,
  record: &RunRecord,
  config: &Map<String, Value>,
  inputs: &[RunInput],
  hook: &ConfigHook,
) -> Result<String, Report> {
  let command = record.config.command();
  let catalog = command_settings(command)?;
  let source = record
    .source_run_id
    .as_ref()
    .filter(|source| manager.store().contains(source))
    .map(|source| manager.store().read(source))
    .transpose()?
    .filter(|source| source.config.command() == command);
  let (baseline, baseline_inputs) = if let Some(source) = &source {
    (source.config.settings()?, compared_inputs(&source.inputs))
  } else {
    let mut defaults = Value::Object(command.config_over_defaults(&Value::Object(Map::new()))?);
    hook(&mut defaults)?;
    let Value::Object(defaults) = defaults else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    (defaults, vec![])
  };
  let differences = setting_differences(
    &catalog,
    &ComparedConfig {
      settings: &baseline,
      inputs: &baseline_inputs,
    },
    &ComparedConfig {
      settings: config,
      inputs: &compared_inputs(inputs),
    },
  );
  let dataset = dataset_label(record.example.as_deref(), inputs);
  Ok(run_name(dataset.as_deref(), &catalog, &differences))
}

fn input_digests(inputs: &[RunInput]) -> Vec<&str> {
  inputs
    .iter()
    .map(|input| input.sha256.as_str())
    .sorted()
    .dedup()
    .collect_vec()
}

fn requested_outputs(config: &Value) -> Option<&Value> {
  config
    .get(OUTPUT_SELECTION_KEY)
    .filter(|selection| !selection.is_null() && selection.as_array().is_none_or(|list| !list.is_empty()))
}
