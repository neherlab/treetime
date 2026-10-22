use crate::command::AppCommandExt;
use crate::command::config_from_settings;
use crate::command::run_command;
use crate::command::{CommandOutcome, OutputFile};
use crate::job::{JobEvent, JobId, JobProgress, JobStarted, TerminalEvent, run_job, run_job_step};
use crate::runs::errors::invalid;
use crate::runs::headline::{RunHeadline, run_headline};
use crate::runs::inputs::hash_inputs;
use crate::runs::manager::{ActiveRun, ConfigHook, RunManager};
use crate::runs::record::{RunError, RunInput, RunRecord, RunStatus};
use crate::runs::warnings::WarningCollector;
use app_commands::command::AppCommand;
use app_commands::json_value::SparseConfig;
use chrono::Utc;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use log::{error, info};
use serde_json::Value;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;
use treetime::cancel::Cancel;
use treetime::progress::RunWarning;
use treetime_utils::io::json::{JsonPretty, json_write_str};

const NO_ERROR_LINES: &str = "the run worker stopped without an error message";

pub fn compute_run(job: &RunJob, cancel: &dyn Cancel, emit: impl Fn(JobEvent) + Send + Sync) -> RunOutcome {
  let progress = JobProgress::new(emit);
  let log = WarningCollector::new(&progress);
  let terminal = run_job(&job.id, cancel, || {
    let prepared = job
      .command
      .prepare_run(&Value::Object(job.config.0.clone()), &job.out_dir)?;
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
        let outcome = compute_run(&prepared.job, &prepared.run.token, |event| prepared.emit(event));
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
    let out_dir = manager.store().out_dir(&id);
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
      let mut config = Value::Object(record.config.settings()?);
      hook(&mut config)?;
      let prepared = command.prepare_run(&config, &out_dir)?;
      let hashed = hash_inputs(command, &prepared.config)?;
      let resolved = config_from_settings(command, &Value::Object(prepared.config.clone()))?;
      let inputs = job_inputs(&hashed.inputs);
      manager.modify(&id, |record| {
        record.config = resolved;
        record.changed_settings = prepared.changed_settings.clone();
        record.inputs = hashed.inputs;
        record.config_hash = Some(hashed.config_hash);
        Ok(())
      })?;
      run.token.check()?;
      Ok(RunJob {
        id: id.clone(),
        command,
        config: SparseConfig(prepared.config),
        out_dir: out_dir.clone(),
        inputs,
      })
    });
    match job {
      Ok(job) => Prepared::Ready(PreparedRun {
        manager,
        run,
        clock,
        job,
      }),
      Err(terminal) => {
        let ended = RunEnd {
          manager: &manager,
          run: &run,
          id: &id,
          out_dir: &out_dir,
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
  Ready(PreparedRun),
  Ended(TerminalEvent),
}

pub struct PreparedRun {
  manager: Arc<RunManager>,
  run: Arc<ActiveRun>,
  clock: Instant,
  job: RunJob,
}

impl PreparedRun {
  pub fn job(&self) -> &RunJob {
    &self.job
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
      out_dir: &self.job.out_dir,
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
  pub out_dir: PathBuf,
  pub inputs: Vec<JobInput>,
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

struct RunEnd<'a> {
  manager: &'a RunManager,
  run: &'a ActiveRun,
  id: &'a JobId,
  out_dir: &'a Path,
  clock: Instant,
}

impl RunEnd<'_> {
  fn record(&self, outcome: RunOutcome) -> TerminalEvent {
    let RunOutcome { terminal, warnings } = outcome;
    let duration = self.clock.elapsed().as_secs_f64();
    let finished = self.manager.modify(self.id, |record| {
      record.finished_at = Some(Utc::now());
      record.warnings = warnings;
      record.duration_seconds = Some(duration);
      match &terminal {
        TerminalEvent::Ok { result, .. } => {
          record.status = RunStatus::Ok;
          record.output_files = relative_outputs(result, self.out_dir);
          record.headline = run_headline(record, self.out_dir).unwrap_or_else(|err| {
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
      error!("When recording the end of run {}: {err:#}", self.id.as_str());
    }
    append_event(self.run, self.id, JobEvent::Terminal { data: terminal.clone() });
    self.manager.finish(self.id);
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

fn relative_outputs(outcome: &CommandOutcome, out_dir: &Path) -> Vec<OutputFile> {
  outcome
    .output_files
    .iter()
    .map(|file| OutputFile {
      path: file.path.strip_prefix(out_dir).unwrap_or(&file.path).to_path_buf(),
      ..file.clone()
    })
    .collect_vec()
}
