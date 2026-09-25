use crate::runs::errors::{InvalidRunRequest, RunConflict, RunNotFound, UploadTooLarge};
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::any::Any;

/// Error of a back-end operation, as the web server and the desktop back end report it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct ErrorResponse {
  /// Class of the error.
  pub code: ErrorCode,
  /// What failed.
  pub message: String,
  /// Underlying causes, outermost first.
  pub causes: Vec<String>,
}

impl ErrorResponse {
  pub fn from_report(report: &Report) -> Self {
    let mut chain = report.chain().map(ToString::to_string);
    Self {
      code: ErrorCode::of(report),
      message: chain.next().unwrap_or_default(),
      causes: chain.collect(),
    }
  }

  pub fn from_panic(payload: &(dyn Any + Send)) -> Self {
    let detail = payload
      .downcast_ref::<&str>()
      .map(|message| (*message).to_owned())
      .or_else(|| payload.downcast_ref::<String>().cloned());
    Self {
      code: ErrorCode::InternalError,
      message: "the back end stopped the operation after an internal error".to_owned(),
      causes: detail.into_iter().collect(),
    }
  }
}

/// Class of a back-end error. The web server answers each class with its own HTTP status.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum ErrorCode {
  /// The run or file does not exist.
  NotFound,
  /// The inputs of a run exceed the upload limit.
  UploadTooLarge,
  /// The run is in a state that does not allow the operation.
  Conflict,
  /// The request is malformed or names an invalid value.
  InvalidRequest,
  /// The back end failed.
  InternalError,
}

impl ErrorCode {
  pub fn of(report: &Report) -> Self {
    if report.downcast_ref::<RunNotFound>().is_some() {
      Self::NotFound
    } else if report.downcast_ref::<UploadTooLarge>().is_some() {
      Self::UploadTooLarge
    } else if report.downcast_ref::<RunConflict>().is_some() {
      Self::Conflict
    } else if report.downcast_ref::<InvalidRunRequest>().is_some()
      || report.downcast_ref::<serde_json::Error>().is_some()
    {
      Self::InvalidRequest
    } else {
      Self::InternalError
    }
  }
}
