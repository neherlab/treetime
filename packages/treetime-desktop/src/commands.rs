use napi::Task;
use napi_derive::napi;
use treetime::commands::ancestral::args::TreetimeAncestralArgs;
use treetime::commands::ancestral::run::run_ancestral_reconstruction;
use treetime::commands::clock::args::TreetimeClockArgs;
use treetime::commands::clock::run::run_clock;
use treetime::commands::mugration::args::TreetimeMugrationArgs;
use treetime::commands::mugration::run::run_mugration;
use treetime::commands::optimize::args::TreetimeOptimizeArgs;
use treetime::commands::optimize::run::run_optimize;
use treetime::commands::prune::args::TreetimePruneArgs;
use treetime::commands::prune::run::run_prune;
use treetime::commands::timetree::args::TreetimeTimetreeArgs;
use treetime::commands::timetree::run::run_timetree_estimation;

fn eyre_to_napi(err: eyre::Report) -> napi::Error {
  napi::Error::new(napi::Status::GenericFailure, format!("{err:?}"))
}

fn json_to_napi(err: serde_json::Error) -> napi::Error {
  napi::Error::new(napi::Status::InvalidArg, format!("{err}"))
}

pub struct AncestralTask {
  args: TreetimeAncestralArgs,
}

impl Task for AncestralTask {
  type Output = ();
  type JsValue = ();

  fn compute(&mut self) -> napi::Result<Self::Output> {
    run_ancestral_reconstruction(&self.args).map_err(eyre_to_napi)
  }

  fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
    Ok(())
  }
}

#[napi]
pub fn ancestral(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<AncestralTask>> {
  let args: TreetimeAncestralArgs = serde_json::from_str(&args_json).map_err(json_to_napi)?;
  Ok(napi::bindgen_prelude::AsyncTask::new(AncestralTask { args }))
}

pub struct ClockTask {
  args: TreetimeClockArgs,
}

impl Task for ClockTask {
  type Output = ();
  type JsValue = ();

  fn compute(&mut self) -> napi::Result<Self::Output> {
    run_clock(&self.args).map_err(eyre_to_napi)?;
    Ok(())
  }

  fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
    Ok(())
  }
}

#[napi]
pub fn clock(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<ClockTask>> {
  let args: TreetimeClockArgs = serde_json::from_str(&args_json).map_err(json_to_napi)?;
  Ok(napi::bindgen_prelude::AsyncTask::new(ClockTask { args }))
}

pub struct TimetreeTask {
  args: TreetimeTimetreeArgs,
}

impl Task for TimetreeTask {
  type Output = ();
  type JsValue = ();

  fn compute(&mut self) -> napi::Result<Self::Output> {
    run_timetree_estimation(&self.args).map_err(eyre_to_napi)
  }

  fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
    Ok(())
  }
}

#[napi]
pub fn timetree(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<TimetreeTask>> {
  let args: TreetimeTimetreeArgs = serde_json::from_str(&args_json).map_err(json_to_napi)?;
  Ok(napi::bindgen_prelude::AsyncTask::new(TimetreeTask { args }))
}

pub struct MugrationTask {
  args: TreetimeMugrationArgs,
}

impl Task for MugrationTask {
  type Output = ();
  type JsValue = ();

  fn compute(&mut self) -> napi::Result<Self::Output> {
    run_mugration(&self.args).map_err(eyre_to_napi)
  }

  fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
    Ok(())
  }
}

#[napi]
pub fn mugration(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<MugrationTask>> {
  let args: TreetimeMugrationArgs = serde_json::from_str(&args_json).map_err(json_to_napi)?;
  Ok(napi::bindgen_prelude::AsyncTask::new(MugrationTask { args }))
}

pub struct OptimizeTask {
  args: TreetimeOptimizeArgs,
}

impl Task for OptimizeTask {
  type Output = ();
  type JsValue = ();

  fn compute(&mut self) -> napi::Result<Self::Output> {
    run_optimize(&self.args).map_err(eyre_to_napi)
  }

  fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
    Ok(())
  }
}

#[napi]
pub fn optimize(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<OptimizeTask>> {
  let args: TreetimeOptimizeArgs = serde_json::from_str(&args_json).map_err(json_to_napi)?;
  Ok(napi::bindgen_prelude::AsyncTask::new(OptimizeTask { args }))
}

pub struct PruneTask {
  args: TreetimePruneArgs,
}

impl Task for PruneTask {
  type Output = ();
  type JsValue = ();

  fn compute(&mut self) -> napi::Result<Self::Output> {
    run_prune(&self.args).map_err(eyre_to_napi)
  }

  fn resolve(&mut self, _env: napi::Env, _output: ()) -> napi::Result<()> {
    Ok(())
  }
}

#[napi]
pub fn prune(args_json: String) -> napi::Result<napi::bindgen_prelude::AsyncTask<PruneTask>> {
  let args: TreetimePruneArgs = serde_json::from_str(&args_json).map_err(json_to_napi)?;
  Ok(napi::bindgen_prelude::AsyncTask::new(PruneTask { args }))
}
