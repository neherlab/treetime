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
use crate::commands::shared::output_args::{
  AncestralOutputSelection, ClockOutputSelection, MugrationOutputSelection, OptimizeOutputSelection,
  PruneOutputSelection, TimetreeOutputSelection,
};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
use crate::commands::timetree::run::run_timetree_estimation;
use crate::config::load::{check_command_config, load_config_document};
use crate::config::properties::{PathRole, leaf_properties};
use crate::config::schema::command_schema;
use crate::config::settings::{remove_setting, setting_ref};
use crate::config::source::{ConfigProblem, ConfigSource, InvalidConfig};
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
#[cfg(feature = "clap")]
use clap::{Command, CommandFactory};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::{JsonSchema, Schema};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::fs;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString, IntoStaticStr, VariantNames};
use treetime::cancel::Cancel;
use treetime::progress::ProgressSink;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::make_error;

const CDS_PLACEHOLDER: &str = "{cds}";

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

  #[cfg(feature = "clap")]
  pub fn cli_command(self) -> Command {
    let mut command = match self {
      Self::Timetree => TreetimeTimetreeArgsRaw::command(),
      Self::Optimize => TreetimeOptimizeArgsRaw::command(),
      Self::Prune => TreetimePruneArgsRaw::command(),
      Self::Ancestral => TreetimeAncestralArgsRaw::command(),
      Self::Clock => TreetimeClockArgsRaw::command(),
      Self::Mugration => TreetimeMugrationArgsRaw::command(),
    };
    command.build();
    command
  }

  pub fn prepare_text(self, source_name: &str, text: &str) -> Result<PreparedCommand, Report> {
    self.prepare_source(source_name, text, None)
  }

  pub fn prepare_value(self, config: &Value) -> Result<PreparedCommand, Report> {
    let text = json_write_str(config, JsonPretty(true))?;
    self.prepare_text("config.json", &text)
  }

  pub fn prepare_run(self, config: &Value, out_dir: &Path) -> Result<PreparedCommand, Report> {
    let mut config = config.clone();
    let Value::Object(settings) = &mut config else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    for leaf in leaf_properties(self.config_schema().as_value())? {
      if leaf.path_role == Some(PathRole::Output) {
        remove_setting(settings, &leaf.key_path);
      }
    }
    let text = json_write_str(&config, JsonPretty(true))?;
    self.prepare_source("config.json", &text, Some(out_dir))
  }

  fn prepare_source(self, source_name: &str, text: &str, run_out: Option<&Path>) -> Result<PreparedCommand, Report> {
    let source = ConfigSource::new(source_name, text);
    match self {
      Self::Timetree => prepare::<TreetimeTimetreeArgsRaw>(&source, text, run_out),
      Self::Optimize => prepare::<TreetimeOptimizeArgsRaw>(&source, text, run_out),
      Self::Prune => prepare::<TreetimePruneArgsRaw>(&source, text, run_out),
      Self::Ancestral => prepare::<TreetimeAncestralArgsRaw>(&source, text, run_out),
      Self::Clock => prepare::<TreetimeClockArgsRaw>(&source, text, run_out),
      Self::Mugration => prepare::<TreetimeMugrationArgsRaw>(&source, text, run_out),
    }
  }
}

pub struct PreparedCommand {
  pub config: Map<String, Value>,
  pub changed_settings: Vec<String>,
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
      output_files: written_files(&outputs)?,
    })
  }
}

/// Result of a command that ran to completion.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct CommandOutcome {
  /// Command that ran.
  pub command: AppCommand,
  /// Files of the command's output plan that exist after the run, sorted by path.
  pub output_files: Vec<OutputFile>,
}

/// One file a command wrote.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
pub struct OutputFile {
  /// Path of the file.
  pub path: PathBuf,
  /// Output selection that produced the file.
  pub kind: OutputSelection,
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

  fn set_run_outputs(&mut self, out_dir: &Path);
}

macro_rules! impl_raw_config {
  ($raw:ty, $args:ty, $variant:ident, $kind:ident, [$($required:expr),* $(,)?]) => {
    impl RawConfig for $raw {
      type Args = $args;

      fn wrap(args: Self::Args) -> CommandArgs {
        CommandArgs::$variant(Box::new(args))
      }

      fn set_run_outputs(&mut self, out_dir: &Path) {
        self.output.output_all = Some(out_dir.to_path_buf());
        self.output_selection = run_output_selection(CommandKind::$kind, &self.output_selection, &[$($required),*]);
      }
    }
  };
}

impl_raw_config!(
  TreetimeTimetreeArgsRaw,
  TreetimeTimetreeArgs,
  Timetree,
  Timetree,
  [TimetreeOutputSelection::Auspice, TimetreeOutputSelection::Tracelog,]
);
impl_raw_config!(
  TreetimeOptimizeArgsRaw,
  TreetimeOptimizeArgs,
  Optimize,
  Optimize,
  [OptimizeOutputSelection::Auspice]
);
impl_raw_config!(
  TreetimePruneArgsRaw,
  TreetimePruneArgs,
  Prune,
  Prune,
  [PruneOutputSelection::Auspice]
);
impl_raw_config!(
  TreetimeAncestralArgsRaw,
  TreetimeAncestralArgs,
  Ancestral,
  Ancestral,
  [AncestralOutputSelection::Auspice]
);
impl_raw_config!(
  TreetimeClockArgsRaw,
  TreetimeClockArgs,
  Clock,
  Clock,
  [ClockOutputSelection::Auspice]
);
impl_raw_config!(
  TreetimeMugrationArgsRaw,
  TreetimeMugrationArgs,
  Mugration,
  Mugration,
  [MugrationOutputSelection::Auspice]
);

fn prepare<R: RawConfig>(source: &ConfigSource, text: &str, run_out: Option<&Path>) -> Result<PreparedCommand, Report> {
  let schema = command_schema::<R>();
  let merged = load_config_document::<R>(source, text)?;
  check_command_config(source, &merged, &schema)?;
  let mut raw: R = serde_json::from_value(merged)?;
  let mut defaults = R::default();
  if let Some(out_dir) = run_out {
    raw.set_run_outputs(out_dir);
    defaults.set_run_outputs(out_dir);
  }
  let config = settings_map(&raw)?;
  let changed_settings = changed_settings(&config, &settings_map(&defaults)?, schema.as_value())?;
  let args = R::Args::try_from(raw)?;
  Ok(PreparedCommand {
    config,
    changed_settings,
    args: R::wrap(args),
  })
}

fn settings_map<R: Serialize>(raw: &R) -> Result<Map<String, Value>, Report> {
  match serde_json::to_value(raw)? {
    Value::Object(settings) => Ok(settings),
    _ => make_error!("a command configuration must serialize to a mapping of settings"),
  }
}

fn changed_settings(
  config: &Map<String, Value>,
  defaults: &Map<String, Value>,
  schema: &Value,
) -> Result<Vec<String>, Report> {
  Ok(
    leaf_properties(schema)?
      .into_iter()
      .filter(|leaf| leaf.path_role.is_none())
      .filter(|leaf| setting_ref(config, &leaf.key_path) != setting_ref(defaults, &leaf.key_path))
      .map(|leaf| leaf.key_path.join("."))
      .collect(),
  )
}

fn run_output_selection<S: Copy + PartialEq + Into<OutputSelection> + IntoEnumIterator>(
  command: CommandKind,
  chosen: &[S],
  required: &[S],
) -> Vec<S> {
  let chosen = if chosen.is_empty() {
    let defaults = command.default_outputs();
    S::iter()
      .filter(|selection| defaults.contains(&(*selection).into()))
      .collect_vec()
  } else {
    chosen.to_vec()
  };
  if chosen
    .iter()
    .any(|selection| (*selection).into() == OutputSelection::All)
  {
    return chosen;
  }
  let missing = required
    .iter()
    .filter(|selection| !chosen.contains(selection))
    .copied()
    .collect_vec();
  [chosen, missing].concat()
}

fn written_files(outputs: &ResolvedOutputs) -> Result<Vec<OutputFile>, Report> {
  let mut files = vec![];
  for (kind, paths) in outputs.paths_by_selection() {
    for path in paths {
      files.extend(existing_files(&path)?.into_iter().map(|path| OutputFile { path, kind }));
    }
  }
  Ok(files.into_iter().sorted().dedup_by(|a, b| a.path == b.path).collect())
}

fn existing_files(planned: &Path) -> Result<Vec<PathBuf>, Report> {
  let name = planned
    .file_name()
    .map(|name| name.to_string_lossy().into_owned())
    .unwrap_or_default();
  let Some((prefix, suffix)) = name.split_once(CDS_PLACEHOLDER) else {
    return Ok(if planned.is_file() {
      vec![planned.to_path_buf()]
    } else {
      vec![]
    });
  };
  let dir = planned.parent().unwrap_or_else(|| Path::new("."));
  if !dir.is_dir() {
    return Ok(vec![]);
  }
  let mut files = vec![];
  for entry in fs::read_dir(dir).wrap_err_with(|| format!("When listing the outputs in '{}'", dir.display()))? {
    let path = entry?.path();
    let matches = path.file_name().is_some_and(|file| {
      let file = file.to_string_lossy();
      file.len() > prefix.len() + suffix.len() && file.starts_with(prefix) && file.ends_with(suffix)
    });
    if matches && path.is_file() {
      files.push(path);
    }
  }
  Ok(files)
}
