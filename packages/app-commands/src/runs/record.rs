use crate::command::{AppCommand, OutputFile};
use crate::job::JobId;
use crate::runs::headline::RunHeadline;
use chrono::{DateTime, Utc};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::path::PathBuf;
use strum_macros::Display;

/// Durable record of one command run, stored as `run.json` in the run's folder.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct RunRecord {
  /// Identifier of the run, also the name of its folder.
  pub id: JobId,
  /// Title shown in run lists.
  pub title: String,
  /// Command the run executes.
  pub command: AppCommand,
  /// Configuration of the run. Before the run starts this is the submitted configuration; afterwards it is the full
  /// resolved configuration, with every default filled in and the outputs the run layer adds.
  pub config: Map<String, Value>,
  /// State of the run.
  pub status: RunStatus,
  /// Whether the run is pinned in run lists.
  pub pinned: bool,
  /// Time the run was created.
  #[schemars(with = "String")]
  pub created_at: DateTime<Utc>,
  /// Time the computation started.
  #[schemars(with = "Option<String>")]
  pub started_at: Option<DateTime<Utc>>,
  /// Time the run ended.
  #[schemars(with = "Option<String>")]
  pub finished_at: Option<DateTime<Utc>>,
  /// Duration of the computation, in seconds.
  pub duration_seconds: Option<f64>,
  /// Version of TreeTime that ran the command.
  pub treetime_version: String,
  /// Input files the run read.
  pub inputs: Vec<RunInput>,
  /// SHA-256 of the canonical resolved configuration, with each input path replaced by that input's SHA-256 and output
  /// paths removed. Two runs with equal hashes used the same settings on the same input contents.
  pub config_hash: Option<String>,
  /// Setting keys whose values differ from the command defaults, as dot-separated key paths.
  pub changed_settings: Vec<String>,
  /// Key results of a finished run, for run lists.
  pub headline: RunHeadline,
  /// Files the run wrote, with paths relative to the run's `out/` folder.
  pub output_files: Vec<OutputFile>,
  /// The error of a failed run.
  pub error: Option<RunError>,
}

impl RunRecord {
  pub fn summary(&self) -> RunSummary {
    RunSummary {
      id: self.id.clone(),
      title: self.title.clone(),
      command: self.command,
      status: self.status,
      pinned: self.pinned,
      created_at: self.created_at,
      finished_at: self.finished_at,
      duration_seconds: self.duration_seconds,
      config_hash: self.config_hash.clone(),
      changed_settings: self.changed_settings.clone(),
      headline: self.headline.clone(),
    }
  }
}

/// State of a run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, Display)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum RunStatus {
  /// Created and waiting to be started, for example while its inputs upload.
  Created,
  /// The computation is running.
  Running,
  /// The command ran to completion.
  Ok,
  /// The command failed or its configuration was rejected.
  Error,
  /// The run stopped because cancellation was requested.
  Cancelled,
  /// The run stopped because the process that ran it stopped.
  Interrupted,
}

impl RunStatus {
  pub const fn is_finished(self) -> bool {
    matches!(self, Self::Ok | Self::Error | Self::Cancelled | Self::Interrupted)
  }
}

/// One input file of a run.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct RunInput {
  /// Dot-separated key path of the setting that names the file.
  pub setting: String,
  /// Path of the file, as the run read it.
  pub path: PathBuf,
  /// Size of the file in bytes.
  pub size: usize,
  /// SHA-256 of the file contents, as lowercase hexadecimal.
  pub sha256: String,
}

/// Error of a failed run.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct RunError {
  /// The error, as the CLI prints it.
  pub message: String,
  /// The errors that caused `message`, outermost first.
  pub causes: Vec<String>,
}

/// Entry of a run list.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct RunSummary {
  /// Identifier of the run.
  pub id: JobId,
  /// Title of the run.
  pub title: String,
  /// Command the run executes.
  pub command: AppCommand,
  /// State of the run.
  pub status: RunStatus,
  /// Whether the run is pinned.
  pub pinned: bool,
  /// Time the run was created.
  #[schemars(with = "String")]
  pub created_at: DateTime<Utc>,
  /// Time the run ended.
  #[schemars(with = "Option<String>")]
  pub finished_at: Option<DateTime<Utc>>,
  /// Duration of the computation, in seconds.
  pub duration_seconds: Option<f64>,
  /// Hash that identifies runs with the same settings and input contents.
  pub config_hash: Option<String>,
  /// Setting keys whose values differ from the command defaults.
  pub changed_settings: Vec<String>,
  /// Key results of a finished run.
  pub headline: RunHeadline,
}

/// Runs, newest first, and the number of runs computing now.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct RunList {
  /// Runs, newest first.
  pub runs: Vec<RunSummary>,
  /// Number of runs whose computation is running now. Every run shares one pool of processing threads.
  pub active_runs: usize,
}

/// Request to create a run.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CreateRunRequest {
  /// Command the run executes.
  pub command: AppCommand,
  /// Configuration of the command, in the form `treetime <command> --config` reads.
  pub config: Value,
  /// Title of the run. Defaults to the command name.
  #[serde(default)]
  pub title: Option<String>,
  /// Whether to wait for an explicit start instead of starting at once, for example to upload inputs first.
  #[serde(default)]
  pub defer_start: bool,
}

/// Request to start a created run.
#[derive(Clone, Debug, Default, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StartRunRequest {
  /// Configuration that replaces the one given at creation, for example to point at uploaded inputs.
  #[serde(default)]
  pub config: Option<Value>,
}

/// Answer to a cancellation request.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct CancelRunResponse {
  /// Whether cancellation was requested; the run ends with a `cancelled` terminal event.
  pub cancelled: bool,
}

/// Changes to the presentation of a run.
#[derive(Clone, Debug, Default, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct UpdateRunRequest {
  /// New title.
  #[serde(default)]
  pub title: Option<String>,
  /// New pinned state.
  #[serde(default)]
  pub pinned: Option<bool>,
}
