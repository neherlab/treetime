use crate::command::{AppCommand, CommandArgs, command_of, settings_map};
use crate::commands::ancestral::args::TreetimeAncestralArgsRaw;
use crate::commands::clock::args::TreetimeClockArgsRaw;
use crate::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
use crate::commands::mugration::args::TreetimeMugrationArgsRaw;
use crate::commands::optimize::args::TreetimeOptimizeArgsRaw;
use crate::commands::prune::args::TreetimePruneArgsRaw;
use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
use app_output::output_plan::ResolvedOutputs;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::path::Path;

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

macro_rules! each_config {
  ($self:expr, |$config:ident| $body:expr) => {
    match $self {
      Self::Timetree { config: $config } => $body,
      Self::Optimize { config: $config } => $body,
      Self::Prune { config: $config } => $body,
      Self::Ancestral { config: $config } => $body,
      Self::Homoplasy { config: $config } => $body,
      Self::Clock { config: $config } => $body,
      Self::Mugration { config: $config } => $body,
    }
  };
}

impl CommandConfig {
  pub fn command(&self) -> AppCommand {
    each_config!(self, |config| command_of(&**config))
  }

  pub fn args(&self) -> Result<CommandArgs, Report> {
    each_config!(self, |config| CommandArgs::try_from((**config).clone()))
  }

  pub fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    self.args()?.resolve_outputs()
  }

  pub fn output_all(&self) -> Option<&Path> {
    each_config!(self, |config| config.output.output_all.as_deref())
  }

  pub fn settings(&self) -> Result<Map<String, Value>, Report> {
    each_config!(self, |config| settings_map(config))
  }
}
