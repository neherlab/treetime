use crate::command_config::CommandConfig;
use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
use crate::commands::ancestral::run::run_ancestral_reconstruction;
use crate::commands::clock::args::{TreetimeClockArgs, TreetimeClockArgsRaw};
use crate::commands::clock::run::run_clock;
use crate::commands::homoplasy::args::{TreetimeHomoplasyArgs, TreetimeHomoplasyArgsRaw};
use crate::commands::homoplasy::run::run_homoplasy;
use crate::commands::mugration::args::{TreetimeMugrationArgs, TreetimeMugrationArgsRaw};
use crate::commands::mugration::run::run_mugration;
use crate::commands::optimize::args::{TreetimeOptimizeArgs, TreetimeOptimizeArgsRaw};
use crate::commands::optimize::run::run_optimize;
use crate::commands::prune::args::{TreetimePruneArgs, TreetimePruneArgsRaw};
use crate::commands::prune::run::run_prune;
use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
use crate::commands::timetree::run::run_timetree_estimation;
use crate::config::load::read_settings;
use crate::config::schema::command_schema;
use app_output::output_plan::ResolvedOutputs;
#[cfg(feature = "clap")]
use clap::{Command, CommandFactory};
use deser::de::DeserializeOwned;
use deser::{Deserialize, Serialize};
use eyre::Report;
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde_json::{Map, Value};
use std::path::PathBuf;
use strum_macros::{Display, EnumIter, EnumString, IntoStaticStr, VariantNames};
use treetime::cancel::Cancel;
use treetime::progress::{LogSink, StageSink};
use treetime_utils::io::json::to_json_value;
use treetime_utils::make_error;

pub fn command_of<R: RawConfig>(_config: &R) -> AppCommand {
  R::COMMAND
}

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
  JsonSchema,
  Display,
  EnumString,
  EnumIter,
  IntoStaticStr,
  VariantNames,
  Serialize,
  Deserialize,
)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum AppCommand {
  Timetree,
  Clock,
  Ancestral,
  Homoplasy,
  Mugration,
  Optimize,
  Prune,
}

impl AppCommand {
  pub fn config_schema(self) -> Schema {
    self.visit(SchemaOf)
  }

  pub fn config_subschema(self, generator: &mut SchemaGenerator) -> Schema {
    self.visit(SubschemaOf { generator })
  }

  pub fn config_schema_name(self) -> String {
    self.visit(SchemaNameOf)
  }

  pub fn default_config(self) -> Result<Map<String, Value>, Report> {
    self.visit(DefaultsOf)
  }

  pub fn command_config(self, settings: &Value) -> Result<CommandConfig, Report> {
    self.visit(ConfigOf { settings })
  }

  #[cfg(feature = "clap")]
  pub fn cli_command(self) -> Command {
    let mut command = self.visit(CliOf);
    command.build();
    command
  }

  fn visit<V: RawConfigVisitor>(self, visitor: V) -> V::Output {
    match self {
      Self::Timetree => visitor.visit::<TreetimeTimetreeArgsRaw>(),
      Self::Optimize => visitor.visit::<TreetimeOptimizeArgsRaw>(),
      Self::Prune => visitor.visit::<TreetimePruneArgsRaw>(),
      Self::Ancestral => visitor.visit::<TreetimeAncestralArgsRaw>(),
      Self::Homoplasy => visitor.visit::<TreetimeHomoplasyArgsRaw>(),
      Self::Clock => visitor.visit::<TreetimeClockArgsRaw>(),
      Self::Mugration => visitor.visit::<TreetimeMugrationArgsRaw>(),
    }
  }
}

pub enum CommandArgs {
  Timetree(Box<TreetimeTimetreeArgs>),
  Optimize(Box<TreetimeOptimizeArgs>),
  Prune(Box<TreetimePruneArgs>),
  Ancestral(Box<TreetimeAncestralArgs>),
  Homoplasy(Box<TreetimeHomoplasyArgs>),
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
      Self::Homoplasy(_) => AppCommand::Homoplasy,
      Self::Clock(_) => AppCommand::Clock,
      Self::Mugration(_) => AppCommand::Mugration,
    }
  }

  pub fn execute(&self, cancel: &dyn Cancel, stages: &dyn StageSink, log: &dyn LogSink) -> Result<(), Report> {
    match self {
      Self::Timetree(args) => run_timetree_estimation(args, cancel, stages, log),
      Self::Optimize(args) => run_optimize(args, cancel, stages, log),
      Self::Prune(args) => run_prune(args, cancel, stages, log),
      Self::Ancestral(args) => run_ancestral_reconstruction(args, cancel, stages, log),
      Self::Homoplasy(args) => run_homoplasy(args, cancel, stages, log).map(|_| ()),
      Self::Clock(args) => run_clock(args, cancel, stages, log).map(|_| ()),
      Self::Mugration(args) => run_mugration(args, cancel, stages, log),
    }
  }

  pub fn resolve_outputs(&self) -> Result<ResolvedOutputs, Report> {
    match self {
      Self::Timetree(args) => args.outputs.resolve(),
      Self::Optimize(args) => args.outputs.resolve(),
      Self::Prune(args) => args.outputs.resolve(),
      Self::Ancestral(args) => args.outputs.resolve(),
      Self::Homoplasy(args) => args.outputs.resolve(),
      Self::Clock(args) => args.outputs.resolve(),
      Self::Mugration(args) => args.outputs.resolve(),
    }
  }

  pub fn alignments(&self) -> &[PathBuf] {
    match self {
      Self::Timetree(args) => args.alignments(),
      Self::Optimize(args) => args.alignments(),
      Self::Prune(args) => args.alignments(),
      Self::Ancestral(args) => args.alignments(),
      Self::Homoplasy(args) => args.alignments(),
      Self::Clock(_) | Self::Mugration(_) => &[],
    }
  }
}

pub trait RawConfig: Serialize + DeserializeOwned + Default + JsonSchema + Clone + RawCli {
  const COMMAND: AppCommand;

  type Args: TryFrom<Self, Error = Report>;

  fn wrap(args: Self::Args) -> CommandArgs;

  fn into_config(self) -> CommandConfig;

  #[must_use]
  fn with_cli_only_settings_of(self, cli: &Self) -> Self;
}

#[cfg(feature = "clap")]
pub trait RawCli: CommandFactory {}

#[cfg(feature = "clap")]
impl<T: CommandFactory> RawCli for T {}

#[cfg(not(feature = "clap"))]
pub trait RawCli {}

#[cfg(not(feature = "clap"))]
impl<T> RawCli for T {}

trait RawConfigVisitor {
  type Output;

  fn visit<R: RawConfig>(self) -> Self::Output;
}

struct SchemaOf;

impl RawConfigVisitor for SchemaOf {
  type Output = Schema;

  fn visit<R: RawConfig>(self) -> Schema {
    command_schema::<R>()
  }
}

struct SubschemaOf<'a> {
  generator: &'a mut SchemaGenerator,
}

impl RawConfigVisitor for SubschemaOf<'_> {
  type Output = Schema;

  fn visit<R: RawConfig>(self) -> Schema {
    self.generator.subschema_for::<R>()
  }
}

struct SchemaNameOf;

impl RawConfigVisitor for SchemaNameOf {
  type Output = String;

  fn visit<R: RawConfig>(self) -> String {
    R::schema_name().into_owned()
  }
}

struct DefaultsOf;

impl RawConfigVisitor for DefaultsOf {
  type Output = Result<Map<String, Value>, Report>;

  fn visit<R: RawConfig>(self) -> Self::Output {
    settings_map(&R::default())
  }
}

struct ConfigOf<'a> {
  settings: &'a Value,
}

impl RawConfigVisitor for ConfigOf<'_> {
  type Output = Result<CommandConfig, Report>;

  fn visit<R: RawConfig>(self) -> Self::Output {
    Ok(read_settings::<R>(self.settings)?.into_config())
  }
}

#[cfg(feature = "clap")]
struct CliOf;

#[cfg(feature = "clap")]
impl RawConfigVisitor for CliOf {
  type Output = Command;

  fn visit<R: RawConfig>(self) -> Command {
    R::command()
  }
}

macro_rules! impl_raw_config {
  ($raw:ty, $args:ty, $variant:ident) => {
    impl RawConfig for $raw {
      const COMMAND: AppCommand = AppCommand::$variant;

      type Args = $args;

      fn wrap(args: Self::Args) -> CommandArgs {
        CommandArgs::$variant(Box::new(args))
      }

      fn into_config(self) -> CommandConfig {
        CommandConfig::$variant { config: Box::new(self) }
      }

      fn with_cli_only_settings_of(self, cli: &Self) -> Self {
        Self {
          config_args: cli.config_args.clone(),
          jobs_args: cli.jobs_args.clone(),
          ..self
        }
      }
    }

    impl TryFrom<$raw> for CommandArgs {
      type Error = Report;

      fn try_from(raw: $raw) -> Result<Self, Report> {
        Ok(<$raw as RawConfig>::wrap(<$args>::try_from(raw)?))
      }
    }
  };
}

impl_raw_config!(TreetimeTimetreeArgsRaw, TreetimeTimetreeArgs, Timetree);
impl_raw_config!(TreetimeOptimizeArgsRaw, TreetimeOptimizeArgs, Optimize);
impl_raw_config!(TreetimePruneArgsRaw, TreetimePruneArgs, Prune);
impl_raw_config!(TreetimeAncestralArgsRaw, TreetimeAncestralArgs, Ancestral);
impl_raw_config!(TreetimeHomoplasyArgsRaw, TreetimeHomoplasyArgs, Homoplasy);
impl_raw_config!(TreetimeClockArgsRaw, TreetimeClockArgs, Clock);
impl_raw_config!(TreetimeMugrationArgsRaw, TreetimeMugrationArgs, Mugration);

pub fn settings_map<R: Serialize>(raw: &R) -> Result<Map<String, Value>, Report> {
  match to_json_value(&raw)? {
    Value::Object(settings) => Ok(settings),
    _ => make_error!("a command configuration must serialize to a mapping of settings"),
  }
}
