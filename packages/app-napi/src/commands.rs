use crate::progress::{self, NapiProgressSink};
use app_api::datasets::discover_datasets;
use app_api::progress::CancelledError;
use app_api::version::version_info;
use app_api::{
  TreetimeAncestralArgs, TreetimeClockArgs, TreetimeMugrationArgs, TreetimeOptimizeArgs, TreetimePruneArgs,
  TreetimeTimetreeArgs,
};
use napi::threadsafe_function::ThreadsafeFunction;
use napi::Task;
use napi_derive::napi;
use std::path::Path;
use std::sync::Arc;

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
  ($task_name:ident, $args_type:ty, $api_fn:path, $napi_fn:ident) => {
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
        let result = $api_fn(&self.args, &sink).map_err(|e| eyre_to_napi(&e))?;
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

define_task!(
  AncestralTask,
  TreetimeAncestralArgs,
  app_api::commands::ancestral,
  ancestral
);
define_task!(ClockTask, TreetimeClockArgs, app_api::commands::clock, clock);
define_task!(
  TimetreeTask,
  TreetimeTimetreeArgs,
  app_api::commands::timetree,
  timetree
);
define_task!(
  MugrationTask,
  TreetimeMugrationArgs,
  app_api::commands::mugration,
  mugration
);
define_task!(
  OptimizeTask,
  TreetimeOptimizeArgs,
  app_api::commands::optimize,
  optimize
);
define_task!(PruneTask, TreetimePruneArgs, app_api::commands::prune, prune);
