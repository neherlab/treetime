use crate::command::{AppCommand, CommandArgs};
use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
use crate::commands::clock::args::TreetimeClockArgsRaw;
use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
use crate::commands::prune::args::TreetimePruneArgsRaw;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use crate::config::load::check_command_config;
use crate::config::source::ConfigSource;
use crate::json_value::SparseConfig;
use crate::runs::errors::invalid;
use app_output::output_plan::ResolvedOutputs;
use eyre::{Report, WrapErr};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::path::Path;
use treetime_utils::error::report_to_string;
use treetime_utils::make_error;

const REQUEST_SOURCE: &str = "config.json";

/// A command with its complete configuration: every setting, defaults included.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "command", content = "config", rename_all = "kebab-case")]
pub enum CommandConfig {
  Timetree(Box<TreetimeTimetreeArgsRaw>),
  Optimize(Box<TreetimeOptimizeArgsRaw>),
  Prune(Box<TreetimePruneArgsRaw>),
  Ancestral(Box<TreetimeAncestralArgsRaw>),
  Clock(Box<TreetimeClockArgsRaw>),
  Mugration(Box<TreetimeMugrationArgsRaw>),
}

impl CommandConfig {
  pub fn from_request(command: AppCommand, config: &SparseConfig) -> Result<Self, Report> {
    let value = Value::Object(config.0.clone());
    let source = ConfigSource::new(REQUEST_SOURCE, serde_json::to_string_pretty(&value)?);
    check_command_config(&source, &value, &command.config_schema())
      .map_err(|report| invalid(report_to_string(&report)))?;
    Self::from_settings(command, &value).map_err(|report| invalid(report_to_string(&report)))
  }

  pub fn from_settings(command: AppCommand, settings: &Value) -> Result<Self, Report> {
    let settings = Value::Object(command.config_over_defaults(settings)?);
    command
      .command_config(&settings)
      .wrap_err_with(|| format!("When reading the configuration of a `{command}` run"))
  }

  pub fn command(&self) -> AppCommand {
    match self {
      Self::Timetree(_) => AppCommand::Timetree,
      Self::Optimize(_) => AppCommand::Optimize,
      Self::Prune(_) => AppCommand::Prune,
      Self::Ancestral(_) => AppCommand::Ancestral,
      Self::Clock(_) => AppCommand::Clock,
      Self::Mugration(_) => AppCommand::Mugration,
    }
  }

  pub fn args(&self) -> Result<CommandArgs, Report> {
    match self {
      Self::Timetree(config) => CommandArgs::try_from((**config).clone()),
      Self::Optimize(config) => CommandArgs::try_from((**config).clone()),
      Self::Prune(config) => CommandArgs::try_from((**config).clone()),
      Self::Ancestral(config) => CommandArgs::try_from((**config).clone()),
      Self::Clock(config) => CommandArgs::try_from((**config).clone()),
      Self::Mugration(config) => CommandArgs::try_from((**config).clone()),
    }
  }

  pub fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    match self {
      Self::Timetree(config) => config.resolve_outputs(),
      Self::Optimize(config) => config.resolve_outputs(),
      Self::Prune(config) => config.resolve_outputs(),
      Self::Ancestral(config) => config.resolve_outputs(),
      Self::Clock(config) => config.resolve_outputs(),
      Self::Mugration(config) => config.resolve_outputs(),
    }
  }

  pub fn output_all(&self) -> Option<&Path> {
    let output = match self {
      Self::Timetree(config) => &config.output,
      Self::Optimize(config) => &config.output,
      Self::Prune(config) => &config.output,
      Self::Ancestral(config) => &config.output,
      Self::Clock(config) => &config.output,
      Self::Mugration(config) => &config.output,
    };
    output.output_all.as_deref()
  }

  pub fn settings(&self) -> Result<Map<String, Value>, Report> {
    let value = match self {
      Self::Timetree(config) => serde_json::to_value(config),
      Self::Optimize(config) => serde_json::to_value(config),
      Self::Prune(config) => serde_json::to_value(config),
      Self::Ancestral(config) => serde_json::to_value(config),
      Self::Clock(config) => serde_json::to_value(config),
      Self::Mugration(config) => serde_json::to_value(config),
    }?;
    match value {
      Value::Object(settings) => Ok(settings),
      _ => make_error!("a command configuration must serialize to a mapping of settings"),
    }
  }
}
