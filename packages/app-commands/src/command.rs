use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
use crate::commands::ancestral::run::run_ancestral_reconstruction;
use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
use crate::commands::clock::run::run_clock;
use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
use crate::commands::mugration::run::run_mugration;
use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
use crate::commands::optimize::run::run_optimize;
use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
use crate::commands::prune::run::run_prune;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
use crate::commands::timetree::run::run_timetree_estimation;
use crate::config::load::{check_command_config, load_config_document};
use crate::config::schema::command_schema;
use crate::config::source::{ConfigProblem, ConfigSource, InvalidConfig};
use app_output::output_plan::ResolvedOutputs;
use eyre::Report;
use itertools::{Itertools, chain};
use schemars::{JsonSchema, Schema};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::path::PathBuf;
use strum_macros::{Display, EnumIter, EnumString, IntoStaticStr, VariantNames};
use treetime::cancel::Cancel;
use treetime::progress::ProgressSink;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::make_error;

/// Analysis command that every client can run.
#[derive(
  Clone,
  Copy,
  Debug,
  PartialEq,
  Eq,
  PartialOrd,
  Ord,
  Hash,
  Serialize,
  Deserialize,
  JsonSchema,
  Display,
  EnumString,
  EnumIter,
  IntoStaticStr,
  VariantNames,
)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum AppCommand {
  Timetree,
  Optimize,
  Prune,
  Ancestral,
  Clock,
  Mugration,
}

impl AppCommand {
  pub fn config_schema(self) -> Schema {
    match self {
      Self::Timetree => command_schema::<TreetimeTimetreeArgsRaw>(),
      Self::Optimize => command_schema::<TreetimeOptimizeArgsRaw>(),
      Self::Prune => command_schema::<TreetimePruneArgsRaw>(),
      Self::Ancestral => command_schema::<TreetimeAncestralArgsRaw>(),
      Self::Clock => command_schema::<TreetimeClockArgsRaw>(),
      Self::Mugration => command_schema::<TreetimeMugrationArgsRaw>(),
    }
  }

  pub fn prepare_text(self, source_name: &str, text: &str) -> Result<PreparedCommand, Report> {
    let source = ConfigSource::new(source_name, text);
    match self {
      Self::Timetree => prepare::<TreetimeTimetreeArgsRaw>(&source, text),
      Self::Optimize => prepare::<TreetimeOptimizeArgsRaw>(&source, text),
      Self::Prune => prepare::<TreetimePruneArgsRaw>(&source, text),
      Self::Ancestral => prepare::<TreetimeAncestralArgsRaw>(&source, text),
      Self::Clock => prepare::<TreetimeClockArgsRaw>(&source, text),
      Self::Mugration => prepare::<TreetimeMugrationArgsRaw>(&source, text),
    }
  }

  pub fn prepare_value(self, config: &Value) -> Result<PreparedCommand, Report> {
    let text = json_write_str(config, JsonPretty(true))?;
    self.prepare_text("config.json", &text)
  }
}

pub struct PreparedCommand {
  pub config: Map<String, Value>,
  pub args: CommandArgs,
}

pub enum CommandArgs {
  Timetree(Box<TreetimeTimetreeArgs>),
  Optimize(Box<TreetimeOptimizeArgs>),
  Prune(Box<TreetimePruneArgs>),
  Ancestral(Box<TreetimeAncestralArgs>),
  Clock(Box<TreetimeClockArgs>),
  Mugration(Box<TreetimeMugrationArgs>),
}

impl CommandArgs {
  pub const fn command(&self) -> AppCommand {
    match self {
      Self::Timetree(_) => AppCommand::Timetree,
      Self::Optimize(_) => AppCommand::Optimize,
      Self::Prune(_) => AppCommand::Prune,
      Self::Ancestral(_) => AppCommand::Ancestral,
      Self::Clock(_) => AppCommand::Clock,
      Self::Mugration(_) => AppCommand::Mugration,
    }
  }

  pub fn run(&self, cancel: &dyn Cancel, progress: &dyn ProgressSink) -> Result<CommandOutcome, Report> {
    let outputs = match self {
      Self::Timetree(args) => {
        run_timetree_estimation(args, cancel, progress)?;
        args.resolve_outputs()?
      },
      Self::Optimize(args) => {
        run_optimize(args, cancel, progress)?;
        args.resolve_outputs()?
      },
      Self::Prune(args) => {
        run_prune(args, cancel, progress)?;
        args.resolve_outputs()?
      },
      Self::Ancestral(args) => {
        run_ancestral_reconstruction(args, cancel, progress)?;
        args.resolve_outputs()?
      },
      Self::Clock(args) => {
        run_clock(args, cancel, progress)?;
        args.resolve_outputs()?
      },
      Self::Mugration(args) => {
        run_mugration(args, cancel, progress)?;
        args.resolve_outputs()?
      },
    };
    Ok(CommandOutcome {
      command: self.command(),
      output_files: written_files(&outputs),
    })
  }
}

/// Result of a command that ran to completion.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct CommandOutcome {
  /// Command that ran.
  pub command: AppCommand,
  /// Files of the command's output plan that exist after the run, sorted by path.
  pub output_files: Vec<PathBuf>,
}

/// Outcome of checking a configuration without running it.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "status", rename_all = "kebab-case")]
pub enum CheckConfigResponse {
  /// The configuration is accepted; `config` holds it with every default filled in.
  Valid { config: Map<String, Value> },
  /// The configuration is rejected.
  Invalid {
    /// The error, as the CLI prints it.
    message: String,
    /// The errors that caused `message`, outermost first.
    causes: Vec<String>,
    /// Problems found by parsing and by the schema check, empty for other errors.
    problems: Vec<ConfigProblem>,
    /// The problems drawn against the configuration text, when the text could be parsed.
    rendered: Option<String>,
  },
}

/// Request to check a configuration.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckConfigRequest {
  /// Command the configuration is for.
  pub command: AppCommand,
  /// Configuration as YAML or JSON text.
  pub text: String,
}

pub fn check_config(request: &CheckConfigRequest) -> CheckConfigResponse {
  match request.command.prepare_text("config.yaml", &request.text) {
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

trait RawConfig: Serialize + DeserializeOwned + Default + JsonSchema + Clone {
  type Args: TryFrom<Self, Error = Report>;

  fn wrap(args: Self::Args) -> CommandArgs;
}

macro_rules! impl_raw_config {
  ($raw:ty, $args:ty, $variant:ident) => {
    impl RawConfig for $raw {
      type Args = $args;

      fn wrap(args: Self::Args) -> CommandArgs {
        CommandArgs::$variant(Box::new(args))
      }
    }
  };
}

impl_raw_config!(TreetimeTimetreeArgsRaw, TreetimeTimetreeArgs, Timetree);
impl_raw_config!(TreetimeOptimizeArgsRaw, TreetimeOptimizeArgs, Optimize);
impl_raw_config!(TreetimePruneArgsRaw, TreetimePruneArgs, Prune);
impl_raw_config!(TreetimeAncestralArgsRaw, TreetimeAncestralArgs, Ancestral);
impl_raw_config!(TreetimeClockArgsRaw, TreetimeClockArgs, Clock);
impl_raw_config!(TreetimeMugrationArgsRaw, TreetimeMugrationArgs, Mugration);

fn prepare<R: RawConfig>(source: &ConfigSource, text: &str) -> Result<PreparedCommand, Report> {
  let merged = load_config_document::<R>(source, text)?;
  check_command_config(source, &merged, &command_schema::<R>())?;
  let raw: R = serde_json::from_value(merged)?;
  let Value::Object(config) = serde_json::to_value(&raw)? else {
    return make_error!("a command configuration must serialize to a mapping of settings");
  };
  let args = R::Args::try_from(raw)?;
  Ok(PreparedCommand {
    config,
    args: R::wrap(args),
  })
}

fn written_files(outputs: &ResolvedOutputs) -> Vec<PathBuf> {
  chain!(outputs.tree_outputs.values(), outputs.non_tree_outputs.values())
    .filter(|path| path.is_file())
    .cloned()
    .sorted()
    .dedup()
    .collect()
}
