use crate::command::AppCommand;
use crate::config::code::{ConfigCode, config_code};
use crate::config::source::{ConfigProblem, InvalidConfig};
use crate::runs::inputs::hash_inputs;
use crate::runs::manager::ConfigHook;
use eyre::Report;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::path::Path;
use treetime_utils::error::ReportChain;
use treetime_utils::make_error;

pub const RUN_CONFIG_OUTPUT_DIR: &str = "out";

pub fn run_config(request: &RunConfigRequest, confine: ConfigHook) -> RunConfigResponse {
  let resolved = request
    .command
    .prepare_run(&request.config, Path::new(RUN_CONFIG_OUTPUT_DIR))
    .and_then(|prepared| Ok((config_code(request.command, &prepared.config)?, prepared.config)));
  match resolved {
    Ok((code, config)) => {
      let (config_hash, config_hash_error) = match config_hash(request.command, &config, confine) {
        Ok(hash) => (Some(hash), None),
        Err(report) => (None, Some(format!("{report:#}"))),
      };
      RunConfigResponse::Valid {
        config,
        code,
        config_hash,
        config_hash_error,
      }
    },
    Err(report) => {
      let ReportChain { message, causes } = ReportChain::of(&report);
      RunConfigResponse::Invalid {
        message,
        causes,
        problems: InvalidConfig::problems_of(&report),
      }
    },
  }
}

/// Request to resolve a configuration as a run resolves it, without running it.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RunConfigRequest {
  /// Command the configuration is for.
  pub command: AppCommand,
  /// Configuration of the command, in the form `treetime <command> --config` reads.
  pub config: Value,
}

/// Outcome of resolving a configuration as a run resolves it.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "status", rename_all = "kebab-case")]
pub enum RunConfigResponse {
  /// The configuration is accepted.
  Valid {
    /// The configuration as the run records it, with every default filled in, the outputs the run layer adds, and
    /// `output_all` set to `out`.
    config: Map<String, Value>,
    /// The command line and the YAML config that reproduce the run.
    code: ConfigCode,
    /// Hash that a run of this configuration records, to find finished runs with the same settings and input
    /// contents; absent when an input cannot be read.
    config_hash: Option<String>,
    /// Why `config_hash` is absent.
    config_hash_error: Option<String>,
  },
  /// The configuration is rejected.
  Invalid {
    /// The error, as the CLI prints it.
    message: String,
    /// The errors that caused `message`, outermost first.
    causes: Vec<String>,
    /// Problems found by parsing and by the schema check, empty for other errors.
    problems: Vec<ConfigProblem>,
  },
}

fn config_hash(command: AppCommand, config: &Map<String, Value>, confine: ConfigHook) -> Result<String, Report> {
  let mut confined = Value::Object(config.clone());
  confine(&mut confined)?;
  let Value::Object(confined) = confined else {
    return make_error!("a command configuration must be a mapping of settings");
  };
  Ok(hash_inputs(command, &confined)?.config_hash)
}
