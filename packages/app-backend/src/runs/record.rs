use crate::command::OutputFile;
use crate::job::JobId;
use crate::runs::headline::RunHeadline;
use crate::runs::input_store::StoredInput;
use crate::runs::run_paths::INPUTS_DIR;
use app_commands::command::AppCommand;
use app_commands::command_config::CommandConfig;
use app_commands::json_value::SparseConfig;
use chrono::{DateTime, Utc};
use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use std::path::{Component, PathBuf};
use strum_macros::Display;
use treetime::progress::RunWarning;
use treetime_schema::skip_serializing_optionals;

/// Durable record of one command run, stored as `run.json` in the run's folder. Every path in the record that lies in
/// the run's folder or in the store of uploaded inputs is relative to the run's folder: uploads appear as
/// `inputs/<name>`, or `inputs/<n>/<name>` when two different uploads share a name, and outputs below `out/`. Inputs
/// read from elsewhere keep their absolute paths.
#[derive(Clone, Debug, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct RunRecord {
  /// Identifier of the run, also the name of its folder.
  pub id: JobId,
  /// Title shown in run lists: `Run <local date and time of creation>` until the user renames the run.
  pub title: String,
  /// Command of the run and its configuration, with every default filled in. After the run starts, the configuration
  /// also holds the outputs the run layer adds and the seed the run used.
  #[schemars(flatten)]
  #[deser(flatten)]
  pub config: CommandConfig,
  /// State of the run.
  pub status: RunStatus,
  /// Whether the run is pinned in run lists.
  pub pinned: bool,
  /// Run whose settings this run started from, such as the source of a variant. The source may have been deleted
  /// since.
  pub source_run_id: Option<JobId>,
  /// Group of related runs in the run history: the group of the source run, or for a run without a source, the group
  /// of the earliest run with the same input contents, else the run's own id.
  pub group_id: JobId,
  /// Size of the run in bytes: its output files and the uploaded inputs that no other run references, as of the
  /// last read of the record.
  pub size: usize,
  /// Seed of the random number generator that the run used, for a run with a random step: the configured seed, or
  /// the seed the run drew when the configuration sets none.
  pub seed: Option<u64>,
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
  /// Input files the run read. Before the run starts, the uploaded inputs its configuration names.
  pub inputs: Vec<RunInput>,
  /// Files uploaded for the run before it started; the run drops them when it starts.
  pub uploads: Vec<StoredInput>,
  /// SHA-256 of the canonical resolved configuration, with each input path replaced by that input's SHA-256 and output
  /// paths removed. Two runs with equal hashes used the same settings on the same input contents.
  pub config_hash: Option<String>,
  /// Setting keys whose values differ from the command defaults, as dot-separated key paths.
  pub changed_settings: Vec<String>,
  /// Key results of a finished run, for run lists.
  pub headline: RunHeadline,
  /// Files the run wrote, with paths relative to the run's `out/` folder.
  pub output_files: Vec<OutputFile>,
  /// Warnings of the run, in the order raised. A run that failed keeps the warnings raised before the failure.
  pub warnings: Vec<RunWarning>,
  /// The error of a failed run.
  pub error: Option<RunError>,
}

impl RunRecord {
  pub fn summary(&self) -> RunSummary {
    RunSummary {
      id: self.id.clone(),
      title: self.title.clone(),
      command: self.config.command(),
      status: self.status,
      pinned: self.pinned,
      source_run_id: self.source_run_id.clone(),
      group_id: self.group_id.clone(),
      size: self.size,
      created_at: self.created_at,
      finished_at: self.finished_at,
      duration_seconds: self.duration_seconds,
      config_hash: self.config_hash.clone(),
      changed_settings: self.changed_settings.clone(),
      headline: self.headline.clone(),
    }
  }
}

/// One input file of a run.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct RunInput {
  /// Dot-separated key path of the setting that names the file.
  pub setting: String,
  /// Path of the file: `inputs/<name>` or `inputs/<n>/<name>` relative to the run's folder for an uploaded file, the
  /// absolute path for a file read from elsewhere.
  pub path: PathBuf,
  /// Size of the file in bytes.
  pub size: usize,
  /// SHA-256 of the file contents, as lowercase hexadecimal.
  pub sha256: String,
}

impl RunInput {
  pub fn stored(&self) -> Option<StoredInput> {
    let mut components = self.path.components();
    let is_upload = components.next() == Some(Component::Normal(INPUTS_DIR.as_ref()));
    let name = self.path.file_name()?.to_str()?;
    is_upload.then(|| StoredInput {
      sha256: self.sha256.clone(),
      name: name.to_owned(),
    })
  }
}

/// Error of a failed run.
#[derive(Clone, Debug, JsonSchema, Serialize, Deserialize)]
pub struct RunError {
  /// The error, as the CLI prints it.
  pub message: String,
  /// The errors that caused `message`, outermost first.
  pub causes: Vec<String>,
}

/// Runs, newest first, and the number of runs computing now.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct RunList {
  /// Runs, newest first.
  pub runs: Vec<RunSummary>,
  /// Number of runs whose computation is running now. Every run shares one pool of processing threads.
  pub active_runs: usize,
}

/// Entry of a run list.
#[derive(Clone, Debug, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
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
  /// Run whose settings this run started from. The source may have been deleted since.
  pub source_run_id: Option<JobId>,
  /// Group of related runs in the run history.
  pub group_id: JobId,
  /// Size of the run in bytes: its output files and the uploaded inputs that no other run references.
  pub size: usize,
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

/// State of a run.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Display, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum RunStatus {
  /// Created and waiting to be started, for example while its inputs upload.
  Created,
  /// The computation is running.
  Running,
  /// The command ran to completion.
  Ok,
  /// The command failed: its inputs could not be read, or the computation failed.
  Error,
  /// The run stopped because cancellation was requested.
  Cancelled,
  /// The run stopped because the process that ran it stopped.
  Interrupted,
}

/// Request to create a run.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct CreateRunRequest {
  /// Command the run executes.
  pub command: AppCommand,
  /// Configuration of the command, in the form `treetime <command> --config` reads. A configuration that the command
  /// rejects is answered with an `invalid_request` error.
  pub config: SparseConfig,
  /// Whether to wait for an explicit start instead of starting at once, for example to upload inputs first.
  #[schemars(default)]
  #[deser(default)]
  pub defer_start: bool,
  /// Run whose settings the new run starts from, such as the source of a variant.
  #[schemars(default)]
  #[deser(default)]
  pub source_run_id: Option<JobId>,
}

/// Request to start a created run.
#[derive(Clone, Debug, Default, JsonSchema, Deserialize)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct StartRunRequest {
  /// Command that replaces the one given at creation, for example when the user chose another analysis after
  /// uploading its inputs.
  #[schemars(default)]
  #[deser(default)]
  pub command: Option<AppCommand>,
  /// Configuration that replaces the one given at creation, for example to point at uploaded inputs.
  #[schemars(default)]
  #[deser(default)]
  pub config: Option<SparseConfig>,
  /// Run whose settings the run starts from; replaces the one given at creation.
  #[schemars(default)]
  #[deser(default)]
  pub source_run_id: Option<JobId>,
}

/// Answer to a cancellation request.
#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct CancelRunResponse {
  /// Whether cancellation was requested; the run ends with a `cancelled` terminal event.
  pub cancelled: bool,
}

/// Where to save an output file of a run, or the archive of all its outputs.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct SaveRunRequest {
  /// Path of the file relative to the run's `out/` folder. Unset: a zip archive of the whole `out/` folder.
  #[schemars(default)]
  #[deser(default)]
  pub path: Option<String>,
  /// Absolute path of the file to write.
  pub destination: PathBuf,
}

/// Changes to the presentation of a run.
#[derive(Clone, Debug, Default, JsonSchema, Deserialize)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct UpdateRunRequest {
  /// New title.
  #[schemars(default)]
  #[deser(default)]
  pub title: Option<String>,
  /// New pinned state.
  #[schemars(default)]
  #[deser(default)]
  pub pinned: Option<bool>,
}
