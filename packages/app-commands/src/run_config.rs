use crate::command::{AppCommand, CheckConfigResponse};
use crate::config::source::InvalidConfig;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::path::Path;

pub const RUN_CONFIG_OUTPUT_DIR: &str = "out";

/// Request to resolve a configuration as a run resolves it, without running it.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunConfigRequest {
  /// Command the configuration is for.
  pub command: AppCommand,
  /// Configuration of the command, in the form `treetime <command> --config` reads.
  pub config: Value,
}

pub fn run_config(request: &RunConfigRequest) -> CheckConfigResponse {
  match request
    .command
    .prepare_run(&request.config, Path::new(RUN_CONFIG_OUTPUT_DIR))
  {
    Ok(prepared) => CheckConfigResponse::Valid {
      config: prepared.config,
    },
    Err(report) => {
      let invalid = report.downcast_ref::<InvalidConfig>();
      CheckConfigResponse::Invalid {
        message: report.to_string(),
        causes: report.chain().skip(1).map(ToString::to_string).collect(),
        problems: invalid.map(|invalid| invalid.problems.clone()).unwrap_or_default(),
        rendered: invalid.map(|invalid| invalid.rendered.clone()),
      }
    },
  }
}
