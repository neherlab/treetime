use napi::Task;
use napi_derive::napi;
use treetime_api::progress::NoopProgress;
use treetime_api::{
  TreetimeAncestralArgs, TreetimeClockArgs, TreetimeMugrationArgs, TreetimeOptimizeArgs, TreetimePruneArgs,
  TreetimeTimetreeArgs,
};

fn eyre_to_napi(err: eyre::Report) -> napi::Error {
  napi::Error::new(napi::Status::GenericFailure, format!("{err:?}"))
}

fn json_to_napi(err: serde_json::Error) -> napi::Error {
  napi::Error::new(napi::Status::InvalidArg, format!("{err}"))
}

macro_rules! define_task {
  ($task_name:ident, $args_type:ty, $api_fn:path, $napi_fn:ident) => {
    pub struct $task_name {
      args: $args_type,
    }

    impl Task for $task_name {
      type Output = ();
      type JsValue = ();

      fn compute(&mut self) -> napi::Result<Self::Output> {
        $api_fn(&self.args, &NoopProgress).map_err(eyre_to_napi)
      }

      fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
        Ok(())
      }
    }

    #[napi]
    pub fn $napi_fn(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<$task_name>> {
      let args: $args_type = serde_json::from_str(&args_json).map_err(json_to_napi)?;
      Ok(napi::bindgen_prelude::AsyncTask::new($task_name { args }))
    }
  };
}

define_task!(
  AncestralTask,
  TreetimeAncestralArgs,
  treetime_api::commands::ancestral,
  ancestral
);
define_task!(ClockTask, TreetimeClockArgs, treetime_api::commands::clock, clock);
define_task!(
  TimetreeTask,
  TreetimeTimetreeArgs,
  treetime_api::commands::timetree,
  timetree
);
define_task!(
  MugrationTask,
  TreetimeMugrationArgs,
  treetime_api::commands::mugration,
  mugration
);
define_task!(
  OptimizeTask,
  TreetimeOptimizeArgs,
  treetime_api::commands::optimize,
  optimize
);
define_task!(PruneTask, TreetimePruneArgs, treetime_api::commands::prune, prune);
