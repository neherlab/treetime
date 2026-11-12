use crate::command::AppCommandExt;
use crate::config::catalog::command_settings;
use crate::config::code::{ConfigCode, SettingCommandLine, config_code, setting_command_lines};
use crate::runs::input_store::InputStore;
use crate::runs::inputs::hash_inputs;
use crate::runs::manager::ConfigHook;
use crate::runs::run_paths::PortablePaths;
use app_commands::command::AppCommand;
use app_commands::config::source::{ConfigProblem, InvalidConfig};
use app_commands::json_value::SparseConfig;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::path::Path;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::error::ReportChain;
use treetime_utils::make_error;

pub const RUN_CONFIG_OUTPUT_DIR: &str = "out";

pub fn run_config(request: &RunConfigRequest, confine: &ConfigHook, store: &InputStore) -> RunConfigResponse {
  let command = request.command;
  let resolved = command
    .prepare_run(
      &Value::Object(request.config.0.clone()),
      Path::new(RUN_CONFIG_OUTPUT_DIR),
    )
    .and_then(|prepared| {
      let mut config = prepared.config.clone();
      PortablePaths::of_settings(command, &config, store)?.apply(command, &mut config)?;
      let seed = prepared.args.random_step().and_then(|_| prepared.args.seed());
      let command_lines = setting_command_lines(command, &command_settings(command)?.settings, &config)?;
      Ok((
        config_code(command, &config)?,
        command_lines,
        config,
        prepared.config,
        seed,
      ))
    });
  match resolved {
    Ok((code, command_lines, config, requested, seed)) => {
      let (config_hash, config_hash_error) = match config_hash(command, &requested, confine) {
        Ok(hash) => (Some(hash), None),
        Err(report) => (None, Some(format!("{report:#}"))),
      };
      RunConfigResponse::Valid {
        config: SparseConfig(config),
        code: Box::new(code),
        command_lines,
        config_hash,
        config_hash_error,
        seed,
      }
    },
    Err(report) => RunConfigResponse::invalid(&report),
  }
}

/// Request to resolve a configuration as a run resolves it, without running it.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct RunConfigRequest {
  /// Command the configuration is for.
  pub command: AppCommand,
  /// Configuration of the command, in the form `treetime <command> --config` reads.
  pub config: SparseConfig,
}

/// Outcome of resolving a configuration as a run resolves it.
#[derive(Clone, Debug, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(tag = "status", rename_all = "kebab-case")]
#[deser(tag = "status", rename_all = "kebab-case")]
pub enum RunConfigResponse {
  /// The configuration is accepted.
  Valid {
    /// The configuration as the run records it, with every default filled in, the outputs the run layer adds,
    /// `output_all` set to `out`, and each uploaded input at its path in the run's folder, `inputs/<name>`.
    config: SparseConfig,
    /// The command line and the YAML config that reproduce the run.
    code: Box<ConfigCode>,
    /// The command-line form of each setting of the command for its value in `config`.
    command_lines: Vec<SettingCommandLine>,
    /// Hash that a run of this configuration records, to find finished runs with the same settings and input
    /// contents; absent when an input cannot be read.
    config_hash: Option<String>,
    /// Why `config_hash` is absent.
    config_hash_error: Option<String>,
    /// Seed that the configuration sets for its random step; absent when it sets none, so that the run draws one, and
    /// when the run has no random step. A finished run with the same `config_hash` has the settings of this
    /// configuration when the configuration sets no seed or the seed of that run.
    seed: Option<u64>,
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

impl RunConfigResponse {
  pub fn invalid(report: &Report) -> Self {
    let ReportChain { message, causes } = ReportChain::of(report);
    Self::Invalid {
      message,
      causes,
      problems: InvalidConfig::problems_of(report),
    }
  }
}

fn config_hash(command: AppCommand, config: &Map<String, Value>, confine: &ConfigHook) -> Result<String, Report> {
  let mut confined = Value::Object(config.clone());
  confine(&mut confined)?;
  let Value::Object(confined) = confined else {
    return make_error!("a command configuration must be a mapping of settings");
  };
  Ok(hash_inputs(command, &confined)?.config_hash)
}
