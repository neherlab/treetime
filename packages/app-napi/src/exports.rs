#![allow(
  clippy::expect_used,
  reason = "application layer: counts and indices to f64, integer division for averaging, graph node access and CLI/config setup invariants, and default variant matches"
)]

//! N-API exports.
//!
//! Deserializes each request into its command's openapi-subset request struct, runs the command
//! orchestration in `crate::commands`, and returns the serialized result. Progress and cancellation
//! flow through the process-global threadsafe sinks in `crate::progress`.

use crate::commands::ancestral::{AncestralArgs, run_ancestral};
use crate::commands::clock::{ClockArgs, run_clock};
use crate::commands::mugration::{MugrationArgs, run_mugration};
use crate::commands::optimize::{OptimizeArgs, run_optimize};
use crate::commands::prune::{PruneArgs, run_prune};
use crate::commands::timetree::{TimetreeArgs, run_timetree};
use crate::progress::{self, NapiCancel, NapiProgressSink};
use app_datasets::discover_datasets;
use napi::Task;
use napi::threadsafe_function::ThreadsafeFunction;
use napi_derive::napi;
use std::path::Path;
use std::sync::Arc;
use treetime::cancel::{CancelledError, NoopCancel};
use treetime::progress::NoopProgress;
use treetime_schema::version_info;

#[napi]
pub fn version() -> String {
  serde_json::to_string(&version_info()).expect("version_info serialization failed")
}

#[napi]
pub fn datasets() -> String {
  let data_dir = std::env::var("DATA_DIR").unwrap_or_else(|_| "data".to_owned());
  let datasets = discover_datasets(Path::new(&data_dir));
  serde_json::to_string(&datasets).expect("datasets serialization failed")
}

#[napi]
#[allow(clippy::needless_pass_by_value)]
pub fn ancestral_sync(args_json: String) -> napi::Result<String> {
  let args: AncestralArgs = serde_json::from_str(&args_json).map_err(|e| json_to_napi(&e))?;
  let result = run_ancestral(&args, &NoopCancel, &NoopProgress).map_err(|e| eyre_to_napi(&e))?;
  serde_json::to_string(&result).map_err(|e| json_to_napi(&e))
}

#[napi]
pub fn cancel() {
  progress::cancel();
}

fn eyre_to_napi(err: &eyre::Report) -> napi::Error {
  if err.downcast_ref::<CancelledError>().is_some() {
    napi::Error::new(napi::Status::Cancelled, "Operation cancelled".to_owned())
  } else {
    napi::Error::new(napi::Status::GenericFailure, format!("{err:#}"))
  }
}

fn json_to_napi(err: &serde_json::Error) -> napi::Error {
  napi::Error::new(napi::Status::InvalidArg, format!("{err}"))
}

macro_rules! define_task {
  ($task_name:ident, $args_type:ty, $run_fn:path, $napi_fn:ident) => {
    pub struct $task_name {
      args: $args_type,
      on_event: Arc<ThreadsafeFunction<String, ()>>,
    }

    impl Task for $task_name {
      type Output = String;
      type JsValue = String;

      fn compute(&mut self) -> napi::Result<Self::Output> {
        progress::reset_cancel();
        let sink = NapiProgressSink::new(self.on_event.clone());
        let result = $run_fn(&self.args, &NapiCancel, &sink).map_err(|e| eyre_to_napi(&e))?;
        serde_json::to_string(&result).map_err(|e| json_to_napi(&e))
      }

      fn resolve(&mut self, _env: napi::Env, output: String) -> napi::Result<String> {
        Ok(output)
      }
    }

    #[napi(
      ts_args_type = "argsJson: string, onEvent: (err: Error | null, eventJson: string) => void",
      ts_return_type = "Promise<string>"
    )]
    #[allow(clippy::needless_pass_by_value)]
    pub fn $napi_fn(
      args_json: String,
      on_event: Arc<ThreadsafeFunction<String, ()>>,
    ) -> napi::Result<napi::bindgen_prelude::AsyncTask<$task_name>> {
      let args: $args_type = serde_json::from_str(&args_json).map_err(|e| json_to_napi(&e))?;
      Ok(napi::bindgen_prelude::AsyncTask::new($task_name { args, on_event }))
    }
  };
}

// The async `ancestral` export runs without a threadsafe progress sink: emitting progress through the
// callback triggers a ThreadsafeFunction segfault for this command, so it computes with a no-op sink
// and the callback receives no events. The other commands emit progress normally.
pub struct AncestralTask {
  args: AncestralArgs,
}

impl Task for AncestralTask {
  type Output = String;
  type JsValue = String;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    let result = run_ancestral(&self.args, &NoopCancel, &NoopProgress).map_err(|e| eyre_to_napi(&e))?;
    serde_json::to_string(&result).map_err(|e| json_to_napi(&e))
  }

  fn resolve(&mut self, _env: napi::Env, output: String) -> napi::Result<String> {
    Ok(output)
  }
}

#[napi(
  ts_args_type = "argsJson: string, onEvent: (err: Error | null, eventJson: string) => void",
  ts_return_type = "Promise<string>"
)]
#[allow(clippy::needless_pass_by_value)]
pub fn ancestral(
  args_json: String,
  _on_event: Arc<ThreadsafeFunction<String, ()>>,
) -> napi::Result<napi::bindgen_prelude::AsyncTask<AncestralTask>> {
  let args: AncestralArgs = serde_json::from_str(&args_json).map_err(|e| json_to_napi(&e))?;
  Ok(napi::bindgen_prelude::AsyncTask::new(AncestralTask { args }))
}

define_task!(ClockTask, ClockArgs, run_clock, clock);
define_task!(TimetreeTask, TimetreeArgs, run_timetree, timetree);
define_task!(MugrationTask, MugrationArgs, run_mugration, mugration);
define_task!(OptimizeTask, OptimizeArgs, run_optimize, optimize);
define_task!(PruneTask, PruneArgs, run_prune, prune);
