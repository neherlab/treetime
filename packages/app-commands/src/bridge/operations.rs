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

macro_rules! app_operations {
  ($(
    $name:literal => $variant:ident fn $method:ident($($arg:ident: $ty:ty),*) -> $response:ty;
  )*) => {
    /// Request of an operation of the app back end: the name of the operation and its arguments.
    #[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
    #[serde(tag = "operation", content = "args", deny_unknown_fields)]
    #[allow(
      clippy::large_enum_variant,
      reason = "a request lives for one call; boxing the large configurations would only add allocations"
    )]
    pub enum OperationRequest {
      $(
        #[serde(rename = $name)]
        $variant { $($arg: $ty),* },
      )*
    }

    pub trait Operations {
      $(fn $method(&self, $($arg: $ty),*) -> Result<$response, Report>;)*
    }

    impl OperationRequest {
      pub fn name(&self) -> &'static str {
        match self {
          $(Self::$variant { .. } => $name,)*
        }
      }

      pub fn handle(self, operations: &impl Operations) -> Result<String, Report> {
        match self {
          $(Self::$variant { $($arg),* } => json_write_str(&operations.$method($($arg),*)?, JsonPretty(false)),)*
        }
      }
    }
  };
}

app_operations! {
  "version" => Version fn version() -> VersionInfo;
  "datasets" => Datasets fn datasets() -> DatasetCatalog;
  "check-config" => CheckConfig fn check_config(request: CheckConfigRequest) -> CheckConfigResponse;
  "run-config" => RunConfig fn run_config(request: RunConfigRequest) -> RunConfigResponse;
  "check-inputs" => CheckInputs fn check_inputs(request: CheckInputsRequest) -> InputFacts;
  "list-runs" => ListRuns fn list_runs() -> RunList;
  "create-run" => CreateRun fn create_run(request: CreateRunRequest) -> RunRecord;
  "get-run" => GetRun fn get_run(id: JobId) -> RunRecord;
  "start-run" => StartRun fn start_run(id: JobId, request: StartRunRequest) -> RunRecord;
  "update-run" => UpdateRun fn update_run(id: JobId, request: UpdateRunRequest) -> RunSummary;
  "cancel-run" => CancelRun fn cancel_run(id: JobId) -> CancelRunResponse;
  "delete-run" => DeleteRun fn delete_run(id: JobId) -> ();
  "restore-run" => RestoreRun fn restore_run(id: JobId) -> RunSummary;
  "purge-run" => PurgeRun fn purge_run(id: JobId) -> ();
  "run-files" => RunFiles fn run_files(id: JobId) -> Vec<RunFile>;
  "run-results" => RunResults fn run_results(id: JobId) -> RunResults;
  "run-auspice" => RunAuspice fn run_auspice(id: JobId) -> AuspiceDocument;
  "compare-runs" => CompareRuns fn compare_runs(id: JobId, other: JobId) -> RunComparison;
  "clade-in-runs" => CladeInRuns fn clade_in_runs(request: CladeRequest) -> CladeInRuns;
}
