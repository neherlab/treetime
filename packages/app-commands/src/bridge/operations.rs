use crate::check_config::{CheckConfigRequest, CheckConfigResponse};
use crate::check_inputs::{CheckInputsRequest, InputFacts};
use crate::config::schema::draft2020_generator;
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
use schemars::{JsonSchema, Schema};
use serde::{Deserialize, Serialize};
use strum_macros::IntoStaticStr;
use treetime_schema::VersionInfo;
use treetime_utils::io::json::{JsonPretty, json_write_str};

pub const REQUEST_ARG: &str = "request";

macro_rules! app_operations {
  ($(
    $(#[doc = $doc:literal])*
    $name:literal => $variant:ident fn $method:ident($($arg:ident: $ty:ty),*) -> $response:ty;
    $http:ident $path:literal $operation_id:literal;
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
        $(#[doc = $doc])*
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

    pub fn http_operations() -> Vec<HttpOperation> {
      vec![
        $(HttpOperation {
          name: $name,
          method: HttpMethod::$http,
          path: $path,
          operation_id: $operation_id,
          description: concat!($($doc),*).trim(),
          args: vec![$(OperationArg { name: stringify!($arg), schema: named_schema::<$ty> }),*],
          response: named_schema::<$response>,
        },)*
      ]
    }
  };
}

app_operations! {
  /// Version of TreeTime.
  "version" => Version fn version() -> VersionInfo;
  Get "/api/version" "version";
  /// Example datasets and example configurations in the data directory.
  "datasets" => Datasets fn datasets() -> DatasetCatalog;
  Get "/api/datasets" "datasets";
  /// The configuration with every default filled in, or the problems found in it.
  "check-config" => CheckConfig fn check_config(request: CheckConfigRequest) -> CheckConfigResponse;
  Post "/api/check-config" "configCheck";
  /// The configuration as a run resolves it, with the outputs the run layer adds and the hash the run records, or the
  /// problems found in it.
  "run-config" => RunConfig fn run_config(request: RunConfigRequest) -> RunConfigResponse;
  Post "/api/run-config" "runConfig";
  /// Facts about the input files, read with the readers the commands use.
  "check-inputs" => CheckInputs fn check_inputs(request: CheckInputsRequest) -> InputFacts;
  Post "/api/check-inputs" "inputsCheck";
  /// Runs, newest first, and the number of runs computing now.
  "list-runs" => ListRuns fn list_runs() -> RunList;
  Get "/api/runs" "runsList";
  /// Create a run; it starts at once unless `defer_start` is set.
  "create-run" => CreateRun fn create_run(request: CreateRunRequest) -> RunRecord;
  Post "/api/runs" "runsCreate";
  /// The record of a run.
  "get-run" => GetRun fn get_run(id: JobId) -> RunRecord;
  Get "/api/runs/{id}" "runsGet";
  /// Start a created run.
  "start-run" => StartRun fn start_run(id: JobId, request: StartRunRequest) -> RunRecord;
  Post "/api/runs/{id}/start" "runsStart";
  /// Change the title or pinned state of a run.
  "update-run" => UpdateRun fn update_run(id: JobId, request: UpdateRunRequest) -> RunSummary;
  Patch "/api/runs/{id}" "runsUpdate";
  /// Request cancellation of a run; the run ends with a `cancelled` terminal event.
  "cancel-run" => CancelRun fn cancel_run(id: JobId) -> CancelRunResponse;
  Post "/api/runs/{id}/cancel" "runsCancel";
  /// Move a run to the trash; restore undoes it.
  "delete-run" => DeleteRun fn delete_run(id: JobId) -> ();
  Delete "/api/runs/{id}" "runsDelete";
  /// Bring a run back from the trash.
  "restore-run" => RestoreRun fn restore_run(id: JobId) -> RunSummary;
  Post "/api/runs/{id}/restore" "runsRestore";
  /// Remove a deleted run for good.
  "purge-run" => PurgeRun fn purge_run(id: JobId) -> ();
  Post "/api/runs/{id}/purge" "runsPurge";
  /// Files in the run's `out/` folder, with their sizes and kinds.
  "run-files" => RunFiles fn run_files(id: JobId) -> Vec<RunFile>;
  Get "/api/runs/{id}/files" "runsFiles";
  /// Results of a finished run, read from its output files.
  "run-results" => RunResults fn run_results(id: JobId) -> RunResults;
  Get "/api/runs/{id}/results" "runsResults";
  /// Auspice JSON of a finished run, with the color scales the app displays.
  "run-auspice" => RunAuspice fn run_auspice(id: JobId) -> AuspiceDocument;
  Get "/api/runs/{id}/auspice" "runsAuspice";
  /// Differences of the second run's results from the first run's.
  "compare-runs" => CompareRuns fn compare_runs(id: JobId, other: JobId) -> RunComparison;
  Get "/api/runs/{id}/compare/{other}" "runsCompare";
  /// Nodes of the other finished time-tree runs with the same set of samples below them.
  "clade-in-runs" => CladeInRuns fn clade_in_runs(request: CladeRequest) -> CladeInRuns;
  Post "/api/clade-in-runs" "cladeInRuns";
}

/// An operation as the HTTP server exposes it.
#[derive(Clone, Debug)]
pub struct HttpOperation {
  pub name: &'static str,
  pub method: HttpMethod,
  pub path: &'static str,
  pub operation_id: &'static str,
  pub description: &'static str,
  pub args: Vec<OperationArg>,
  pub response: fn() -> NamedSchema,
}

impl HttpOperation {
  pub fn path_args(&self) -> impl Iterator<Item = &OperationArg> {
    self.args.iter().filter(|arg| arg.name != REQUEST_ARG)
  }

  pub fn body(&self) -> Option<&OperationArg> {
    self.args.iter().find(|arg| arg.name == REQUEST_ARG)
  }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, IntoStaticStr)]
#[strum(serialize_all = "lowercase")]
pub enum HttpMethod {
  Get,
  Post,
  Patch,
  Delete,
}

#[derive(Clone, Debug)]
pub struct OperationArg {
  pub name: &'static str,
  pub schema: fn() -> NamedSchema,
}

pub struct NamedSchema {
  pub name: String,
  pub schema: Schema,
}

impl NamedSchema {
  pub fn is_unit(&self) -> bool {
    self.name == UNIT_SCHEMA_NAME
  }
}

const UNIT_SCHEMA_NAME: &str = "null";

fn named_schema<T: JsonSchema>() -> NamedSchema {
  NamedSchema {
    name: T::schema_name().into_owned(),
    schema: draft2020_generator().into_root_schema_for::<T>(),
  }
}
