use app_commands::job::JobId;
use app_commands::runs::errors::invalid;
use app_commands::runs::manager::{RunManager, StartedRun};
use app_commands::runs::record::{CreateRunRequest, RunRecord};
use eyre::{Report, WrapErr};
use serde_json::Value;
use std::sync::Arc;

pub fn parse_id(id: &str) -> Result<JobId, Report> {
  JobId::parse(id).map_err(|err| invalid(err.to_string()))
}

pub fn create_run(runs: &RunManager, request_json: &str) -> Result<RunRecord, Report> {
  let request: CreateRunRequest = serde_json::from_str(request_json).wrap_err("When reading the run request")?;
  runs.create(request)
}

pub fn start_run(runs: &Arc<RunManager>, id: &str, config_json: Option<&str>) -> Result<StartedRun, Report> {
  let config: Option<Value> = config_json
    .map(serde_json::from_str)
    .transpose()
    .wrap_err("When reading the run configuration")?;
  runs.start(&parse_id(id)?, config, Box::new(|_config: &mut Value| Ok(())))
}
