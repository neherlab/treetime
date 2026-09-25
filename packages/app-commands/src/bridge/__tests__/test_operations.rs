#[cfg(test)]
mod tests {
  use crate::bridge::operations::DesktopRequest;
  use eyre::Report;
  use helpers::FakeBackend;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use serde_json::{Value, json};
  use treetime_utils::assert_error;

  #[test]
  fn test_operations_request_names_the_operation_and_its_arguments() {
    let request: DesktopRequest = serde_json::from_str(r#"{"operation":"get-run","args":{"id":"r1"}}"#).unwrap();
    assert_eq!(
      ("get-run", json!({ "operation": "get-run", "args": { "id": "r1" } })),
      (request.name(), serde_json::to_value(&request).unwrap())
    );
  }

  #[test]
  fn test_operations_handle_answers_with_the_backend_result_as_json() {
    let request: DesktopRequest = serde_json::from_str(r#"{"operation":"cancel-run","args":{"id":"r1"}}"#).unwrap();
    let answer: Value = serde_json::from_str(&request.handle(&FakeBackend).unwrap()).unwrap();
    assert_eq!(json!({ "cancelled": true }), answer);
  }

  #[test]
  fn test_operations_handle_answers_an_operation_without_result_with_null() {
    let request: DesktopRequest = serde_json::from_str(r#"{"operation":"delete-run","args":{"id":"r1"}}"#).unwrap();
    assert_eq!("null", request.handle(&FakeBackend).unwrap());
  }

  #[test]
  fn test_operations_handle_passes_backend_errors_through() {
    let request: DesktopRequest = serde_json::from_str(r#"{"operation":"get-run","args":{"id":"r9"}}"#).unwrap();
    assert_error!(request.handle(&FakeBackend), "no run with id `r9`");
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown_operation(r#"{"operation":"format-disk","args":{}}"#,               "unknown variant `format-disk`, expected one of `version`, `datasets`, `check-config`, `run-config`, `check-inputs`, `list-runs`, `create-run`, `get-run`, `start-run`, `update-run`, `cancel-run`, `delete-run`, `restore-run`, `purge-run`, `run-files`, `run-results`, `compare-runs`, `clade-in-runs` at line 1 column 26")]
  #[case::invalid_id(       r#"{"operation":"get-run","args":{"id":"../x"}}"#,          "invalid job id `../x`: expected 1 to 128 ASCII letters, digits, `-` or `_` at line 1 column 43")]
  #[case::unknown_argument( r#"{"operation":"get-run","args":{"id":"r1","path":"/"}}"#, "unknown field `path`, expected `id` at line 1 column 47")]
  #[case::missing_argument( r#"{"operation":"get-run","args":{}}"#,                     "missing field `id` at line 1 column 32")]
  #[trace]
  fn test_operations_request_rejects_malformed_requests(#[case] json: &str, #[case] expected: &str) {
    assert_error!(serde_json::from_str::<DesktopRequest>(json).map_err(Report::new), expected);
  }

  mod helpers {
    use crate::bridge::operations::DesktopBackend;
    use crate::check_config::{CheckConfigRequest, CheckConfigResponse};
    use crate::check_inputs::{CheckInputsRequest, InputFacts};
    use crate::job::JobId;
    use crate::results::clades::{CladeInRuns, CladeRequest};
    use crate::results::compare::RunComparison;
    use crate::results::run_results::RunResults;
    use crate::run_config::{RunConfigRequest, RunConfigResponse};
    use crate::runs::errors::not_found;
    use crate::runs::files::RunFile;
    use crate::runs::record::{
      CancelRunResponse, CreateRunRequest, RunList, RunRecord, RunSummary, StartRunRequest, UpdateRunRequest,
    };
    use app_datasets::DatasetCatalog;
    use eyre::Report;
    use treetime_schema::{VersionInfo, version_info};
    use treetime_utils::make_error;

    pub(super) struct FakeBackend;

    impl DesktopBackend for FakeBackend {
      fn version(&self) -> Result<VersionInfo, Report> {
        Ok(version_info())
      }

      fn datasets(&self) -> Result<DatasetCatalog, Report> {
        make_error!("datasets is not part of this test")
      }

      fn check_config(&self, _request: CheckConfigRequest) -> Result<CheckConfigResponse, Report> {
        make_error!("check-config is not part of this test")
      }

      fn run_config(&self, _request: RunConfigRequest) -> Result<RunConfigResponse, Report> {
        make_error!("run-config is not part of this test")
      }

      fn check_inputs(&self, _request: CheckInputsRequest) -> Result<InputFacts, Report> {
        make_error!("check-inputs is not part of this test")
      }

      fn list_runs(&self) -> Result<RunList, Report> {
        make_error!("list-runs is not part of this test")
      }

      fn create_run(&self, _request: CreateRunRequest) -> Result<RunRecord, Report> {
        make_error!("create-run is not part of this test")
      }

      fn get_run(&self, id: JobId) -> Result<RunRecord, Report> {
        Err(not_found(format!("no run with id `{}`", id.as_str())))
      }

      fn start_run(&self, _id: JobId, _request: StartRunRequest) -> Result<RunRecord, Report> {
        make_error!("start-run is not part of this test")
      }

      fn update_run(&self, _id: JobId, _request: UpdateRunRequest) -> Result<RunSummary, Report> {
        make_error!("update-run is not part of this test")
      }

      fn cancel_run(&self, _id: JobId) -> Result<CancelRunResponse, Report> {
        Ok(CancelRunResponse { cancelled: true })
      }

      fn delete_run(&self, _id: JobId) -> Result<(), Report> {
        Ok(())
      }

      fn restore_run(&self, _id: JobId) -> Result<RunSummary, Report> {
        make_error!("restore-run is not part of this test")
      }

      fn purge_run(&self, _id: JobId) -> Result<(), Report> {
        make_error!("purge-run is not part of this test")
      }

      fn run_files(&self, _id: JobId) -> Result<Vec<RunFile>, Report> {
        make_error!("run-files is not part of this test")
      }

      fn run_results(&self, _id: JobId) -> Result<RunResults, Report> {
        make_error!("run-results is not part of this test")
      }

      fn compare_runs(&self, _id: JobId, _other: JobId) -> Result<RunComparison, Report> {
        make_error!("compare-runs is not part of this test")
      }

      fn clade_in_runs(&self, _request: CladeRequest) -> Result<CladeInRuns, Report> {
        make_error!("clade-in-runs is not part of this test")
      }
    }
  }
}
