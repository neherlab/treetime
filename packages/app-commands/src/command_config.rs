use crate::command::{AppCommand, CommandArgs};
use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
use crate::commands::clock::args::TreetimeClockArgsRaw;
use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
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
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::path::Path;
use treetime_utils::error::report_to_string;
use treetime_utils::io::json::to_json_value;
use treetime_utils::make_error;

const REQUEST_SOURCE: &str = "config.json";

/// A command with its complete configuration: every setting, defaults included.
#[derive(Clone, Debug, JsonSchema, Serialize, Deserialize)]
#[schemars(tag = "command", rename_all = "kebab-case")]
#[deser(tag = "command", rename_all = "kebab-case")]
pub enum CommandConfig {
  Timetree { config: Box<TreetimeTimetreeArgsRaw> },
  Optimize { config: Box<TreetimeOptimizeArgsRaw> },
  Prune { config: Box<TreetimePruneArgsRaw> },
  Ancestral { config: Box<TreetimeAncestralArgsRaw> },
  Homoplasy { config: Box<TreetimeHomoplasyArgsRaw> },
  Clock { config: Box<TreetimeClockArgsRaw> },
  Mugration { config: Box<TreetimeMugrationArgsRaw> },
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
      Self::Timetree { .. } => AppCommand::Timetree,
      Self::Optimize { .. } => AppCommand::Optimize,
      Self::Prune { .. } => AppCommand::Prune,
      Self::Ancestral { .. } => AppCommand::Ancestral,
      Self::Homoplasy { .. } => AppCommand::Homoplasy,
      Self::Clock { .. } => AppCommand::Clock,
      Self::Mugration { .. } => AppCommand::Mugration,
    }
  }

  pub fn args(&self) -> Result<CommandArgs, Report> {
    match self {
      Self::Timetree { config } => CommandArgs::try_from((**config).clone()),
      Self::Optimize { config } => CommandArgs::try_from((**config).clone()),
      Self::Prune { config } => CommandArgs::try_from((**config).clone()),
      Self::Ancestral { config } => CommandArgs::try_from((**config).clone()),
      Self::Homoplasy { config } => CommandArgs::try_from((**config).clone()),
      Self::Clock { config } => CommandArgs::try_from((**config).clone()),
      Self::Mugration { config } => CommandArgs::try_from((**config).clone()),
    }
  }

  pub fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    match self {
      Self::Timetree { config } => config.resolve_outputs(),
      Self::Optimize { config } => config.resolve_outputs(),
      Self::Prune { config } => config.resolve_outputs(),
      Self::Ancestral { config } => config.resolve_outputs(),
      Self::Homoplasy { config } => config.resolve_outputs(),
      Self::Clock { config } => config.resolve_outputs(),
      Self::Mugration { config } => config.resolve_outputs(),
    }
  }

  pub fn output_all(&self) -> Option<&Path> {
    let output = match self {
      Self::Timetree { config } => &config.output,
      Self::Optimize { config } => &config.output,
      Self::Prune { config } => &config.output,
      Self::Ancestral { config } => &config.output,
      Self::Homoplasy { config } => &config.output,
      Self::Clock { config } => &config.output,
      Self::Mugration { config } => &config.output,
    };
    output.output_all.as_deref()
  }

  pub fn settings(&self) -> Result<Map<String, Value>, Report> {
    let value = match self {
      Self::Timetree { config } => to_json_value(&config),
      Self::Optimize { config } => to_json_value(&config),
      Self::Prune { config } => to_json_value(&config),
      Self::Ancestral { config } => to_json_value(&config),
      Self::Homoplasy { config } => to_json_value(&config),
      Self::Clock { config } => to_json_value(&config),
      Self::Mugration { config } => to_json_value(&config),
    }?;
    match value {
      Value::Object(settings) => Ok(settings),
      _ => make_error!("a command configuration must serialize to a mapping of settings"),
    }
  }
}
