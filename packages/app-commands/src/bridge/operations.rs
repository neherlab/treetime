use crate::check_config::{CheckConfigRequest, CheckConfigResponse};
use crate::check_inputs::{CheckInputsRequest, InputFacts};
use crate::datasets::DatasetCatalog;
use crate::job::JobId;
use crate::results::auspice::AuspiceDocument;
use crate::results::clades::{CladeInRuns, CladeRequest};
use crate::results::compare::RunComparison;
use crate::results::run_results::RunResults;
use crate::run_config::{RunConfigRequest, RunConfigResponse};
use crate::runs::files::RunFile;
use crate::runs::record::{
  CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
};
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use treetime_schema::VersionInfo;
use treetime_utils::io::json::{JsonPretty, json_write_str};

macro_rules! desktop_operations {
  ($(
    $(#[doc = $doc:literal])*
    $name:literal => $variant:ident fn $method:ident($($arg:ident: $ty:ty),*) -> $response:ty;
  )*) => {
    /// Request to the desktop back end: the name of an operation and its arguments.
    #[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
    #[serde(tag = "operation", content = "args", deny_unknown_fields)]
    #[allow(
      clippy::large_enum_variant,
      reason = "a request lives for one call; boxing the large configurations would only add allocations"
    )]
    pub enum DesktopRequest {
      $(
        $(#[doc = $doc])*
        #[serde(rename = $name)]
        $variant { $($arg: $ty),* },
      )*
    }

    pub trait DesktopBackend {
      $(fn $method(&self, $($arg: $ty),*) -> Result<$response, Report>;)*
    }

    impl DesktopRequest {
      pub fn name(&self) -> &'static str {
        match self {
          $(Self::$variant { .. } => $name,)*
        }
      }

      pub fn handle(self, backend: &impl DesktopBackend) -> Result<String, Report> {
        match self {
          $(Self::$variant { $($arg),* } => json_write_str(&backend.$method($($arg),*)?, JsonPretty(false)),)*
        }
      }
    }
  };
}

desktop_operations! {
  /// Version of TreeTime.
  "version" => Version fn version() -> VersionInfo;
  /// Example datasets and example configurations.
  "datasets" => Datasets fn datasets() -> DatasetCatalog;
  /// The configuration with every default filled in, or the problems found in it.
  "check-config" => CheckConfig fn check_config(request: CheckConfigRequest) -> CheckConfigResponse;
  /// The configuration as a run resolves it.
  "run-config" => RunConfig fn run_config(request: RunConfigRequest) -> RunConfigResponse;
  /// Facts about input files, read with the readers the commands use.
  "check-inputs" => CheckInputs fn check_inputs(request: CheckInputsRequest) -> InputFacts;
  /// Runs, newest first.
  "list-runs" => ListRuns fn list_runs() -> RunList;
  /// Create a run; it starts at once unless `defer_start` is set.
  "create-run" => CreateRun fn create_run(request: CreateRunRequest) -> RunRecord;
  /// The record of a run.
  "get-run" => GetRun fn get_run(id: JobId) -> RunRecord;
  /// Start a created run.
  "start-run" => StartRun fn start_run(id: JobId, request: StartRunRequest) -> RunRecord;
  /// Change the title or pinned state of a run.
  "update-run" => UpdateRun fn update_run(id: JobId, request: UpdateRunRequest) -> RunSummary;
  /// Request cancellation of a run.
  "cancel-run" => CancelRun fn cancel_run(id: JobId) -> CancelRunResponse;
  /// Move a run to the trash.
  "delete-run" => DeleteRun fn delete_run(id: JobId) -> ();
  /// Bring a run back from the trash.
  "restore-run" => RestoreRun fn restore_run(id: JobId) -> RunSummary;
  /// Remove a deleted run for good.
  "purge-run" => PurgeRun fn purge_run(id: JobId) -> ();
  /// Files in the run's `out/` folder.
  "run-files" => RunFiles fn run_files(id: JobId) -> Vec<RunFile>;
  /// Results of a finished run, read from its output files.
  "run-results" => RunResults fn run_results(id: JobId) -> RunResults;
  /// Auspice JSON of a finished run, with the color scales the app displays.
  "run-auspice" => RunAuspice fn run_auspice(id: JobId) -> AuspiceDocument;
  /// Differences of the second run's results from the first run's.
  "compare-runs" => CompareRuns fn compare_runs(id: JobId, other: JobId) -> RunComparison;
  /// Nodes of other finished time-tree runs with the same samples below them.
  "clade-in-runs" => CladeInRuns fn clade_in_runs(request: CladeRequest) -> CladeInRuns;
}
