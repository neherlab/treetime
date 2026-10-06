use crate::check_inputs::{CommandInput, InputKind, InputNeed};
use crate::command_config::CommandConfig;
use crate::commands::ancestral::aa_node_data::CDS_PLACEHOLDERS;
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
use crate::commands::shared::output_args::{
  AncestralOutputSelection, ClockOutputSelection, HomoplasyOutputSelection, MugrationOutputSelection,
  OptimizeOutputSelection, PruneOutputSelection, TimetreeOutputSelection,
};
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::timetree::args::{TreetimeTimetreeArgs, TreetimeTimetreeArgsRaw};
use crate::commands::timetree::run::run_timetree_estimation;
#[cfg(feature = "clap")]
use crate::config::cli_rules::cli_rule_diagnostics;
use crate::config::load::{check_command_config, load_config_document, merge_value};
use crate::config::properties::{PathRole, leaf_properties};
use crate::config::schema::command_schema;
use crate::config::settings::{remove_setting, setting_ref};
use crate::config::source::ConfigSource;
#[cfg(feature = "clap")]
use crate::config::source::render_and_bail;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
#[cfg(feature = "clap")]
use clap::{Command, CommandFactory};
use deser::adapters::As;
use deser::de::DeserializeOwned;
use deser::{Deserialize, Serialize};
use deser_serde::Serde;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde_json::{Map, Value};
use std::fs;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString, IntoStaticStr, VariantNames};
use treetime::cancel::Cancel;
use treetime::progress::{LogSink, StageSink};
use treetime_utils::io::json::from_json_value;
use treetime_utils::io::json::to_json_value;
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

  pub fn config_over_defaults(self, settings: &Value) -> Result<Map<String, Value>, Report> {
    let mut config = Value::Object(self.default_config()?);
    if settings.is_object() {
      merge_value(&mut config, settings);
    }
    match config {
      Value::Object(config) => Ok(config),
      _ => make_error!("a command configuration must be a mapping of settings"),
    }
  }

  pub const fn inputs(self) -> &'static [CommandInput] {
    const TREE: CommandInput = CommandInput::new(InputKind::Tree, InputNeed::Required);
    const METADATA: CommandInput = CommandInput::new(InputKind::Metadata, InputNeed::Required);
    const ALIGNMENT: CommandInput = CommandInput::new(InputKind::Alignment, InputNeed::Required);
    const ALIGNMENT_RECOMMENDED: CommandInput = CommandInput::new(InputKind::Alignment, InputNeed::Recommended);
    const ALIGNMENT_OPTIONAL: CommandInput = CommandInput::new(InputKind::Alignment, InputNeed::Optional);
    match self {
      Self::Timetree => &[TREE, METADATA, ALIGNMENT_RECOMMENDED],
      Self::Clock => &[TREE, METADATA, ALIGNMENT_OPTIONAL],
      Self::Ancestral | Self::Homoplasy | Self::Optimize => &[TREE, ALIGNMENT],
      Self::Mugration => &[TREE, METADATA],
      Self::Prune => &[TREE, ALIGNMENT_OPTIONAL],
    }
  }

  pub const fn uses_dates(self) -> bool {
    matches!(self, Self::Timetree | Self::Clock)
  }

  pub const fn title(self) -> &'static str {
    match self {
      Self::Timetree => "Time tree",
      Self::Clock => "Clock signal",
      Self::Ancestral => "Ancestral sequences",
      Self::Homoplasy => "Homoplasy",
      Self::Mugration => "Discrete traits",
      Self::Optimize => "Branch lengths",
      Self::Prune => "Prune tree",
    }
  }

  pub const fn description(self) -> &'static str {
    match self {
      Self::Timetree => "Date the ancestors and estimate the clock rate from sampling dates.",
      Self::Clock => "Root-to-tip regression: clock rate, temporal signal and outliers.",
      Self::Ancestral => "Infer ancestral sequences and the mutations on each branch.",
      Self::Homoplasy => "Mutations that arise on several branches, compared with a Poisson expectation.",
      Self::Mugration => "Ancestral states of a metadata column, such as country or host.",
      Self::Optimize => "Maximum-likelihood branch lengths on a fixed topology.",
      Self::Prune => "Collapse short or empty branches and remove listed samples.",
    }
  }

  pub const fn main_settings(self) -> &'static [&'static str] {
    match self {
      Self::Timetree => &[
        "clock_rate",
        "clock_std_dev",
        "confidence",
        "covariation",
        "coalescent",
        "coalescent_opt",
        "coalescent_skyline",
        "skyline_n_points",
        "skyline_stiffness",
        "reroot",
        "keep_root",
        "clock_filter",
        "relax",
        "keep_polytomies",
        "model",
        "model_params",
        "max_iter",
      ],
      Self::Clock => &[
        "reroot",
        "keep_root",
        "clock_filter",
        "covariation",
        "allow_negative_rate",
        "metadata_id_columns",
      ],
      Self::Ancestral => &[
        "method_anc",
        "model",
        "model_params",
        "gap_fill",
        "reconstruct_tip_states",
      ],
      Self::Homoplasy => &["method_anc", "model", "gap_fill", "constant_sites", "drms", "rescale"],
      Self::Mugration => &["attribute", "weights", "pc", "missing_data", "sampling_bias_correction"],
      Self::Optimize => &["opt_method", "reroot", "divergence_units", "no_indels", "max_iter"],
      Self::Prune => &[
        "prune_short",
        "prune_empty",
        "merge_shared_mutations",
        "prune_nodes_list",
      ],
    }
  }

  #[cfg(feature = "clap")]
  pub fn cli_command(self) -> Command {
    let mut command = self.visit(CliOf);
    command.build();
    command
  }

  pub fn prepare_text(self, source_name: &str, text: &str) -> Result<PreparedCommand, Report> {
    self.prepare_source(source_name, text, None)
  }

  pub fn prepare_value(self, config: &Value) -> Result<PreparedCommand, Report> {
    let text = json_write_str(&As::<_, Serde>::new(config), JsonPretty(true))?;
    self.prepare_text("config.json", &text)
  }

  pub fn prepare_run(self, config: &Value, out_dir: &Path) -> Result<PreparedCommand, Report> {
    let mut config = config.clone();
    let Value::Object(settings) = &mut config else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    self.remove_output_paths(settings)?;
    let text = json_write_str(&As::<_, Serde>::new(&config), JsonPretty(true))?;
    self.prepare_source("config.json", &text, Some(out_dir))
  }

  pub fn remove_output_paths(self, settings: &mut Map<String, Value>) -> Result<(), Report> {
    for leaf in leaf_properties(self.config_schema().as_value())? {
      if leaf.path_role == Some(PathRole::Output) {
        remove_setting(settings, &leaf.key_path);
      }
    }
    Ok(())
  }

  fn prepare_source(self, source_name: &str, text: &str, run_out: Option<&Path>) -> Result<PreparedCommand, Report> {
    let source = ConfigSource::new(source_name, text);
    self.visit(Prepare {
      command: self,
      source: &source,
      text,
      run_out,
    })
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

  pub fn run(&self, cancel: &dyn Cancel, stages: &dyn StageSink, log: &dyn LogSink) -> Result<CommandOutcome, Report> {
    self.execute(cancel, stages, log)?;
    let outputs = match self {
      Self::Timetree(args) => args.resolve_outputs()?,
      Self::Optimize(args) => args.resolve_outputs()?,
      Self::Prune(args) => args.resolve_outputs()?,
      Self::Ancestral(args) => args.resolve_outputs()?,
      Self::Homoplasy(args) => args.resolve_outputs()?,
      Self::Clock(args) => args.resolve_outputs()?,
      Self::Mugration(args) => args.resolve_outputs()?,
    };
    Ok(CommandOutcome {
      command: self.command(),
      output_files: written_files(&outputs)?,
    })
  }
}

/// Result of a command that ran to completion.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct CommandOutcome {
  /// Command that ran.
  pub command: AppCommand,
  /// Files of the command's output plan that exist after the run, sorted by path.
  pub output_files: Vec<OutputFile>,
}

/// One file a command wrote.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, JsonSchema, Serialize, Deserialize)]
pub struct OutputFile {
  /// Path of the file.
  pub path: PathBuf,
  /// Output selection that produced the file.
  pub kind: OutputSelection,
}

trait RawConfig: Serialize + DeserializeOwned + Default + JsonSchema + Clone + RawCli {
  type Args: TryFrom<Self, Error = Report>;

  fn wrap(args: Self::Args) -> CommandArgs;

  fn into_config(self) -> CommandConfig;

  fn set_run_outputs(&mut self, out_dir: &Path);
}

#[cfg(feature = "clap")]
trait RawCli: CommandFactory {}

#[cfg(feature = "clap")]
impl<T: CommandFactory> RawCli for T {}

#[cfg(not(feature = "clap"))]
trait RawCli {}

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
    Ok(from_json_value::<R>(self.settings)?.into_config())
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

struct Prepare<'a> {
  command: AppCommand,
  source: &'a ConfigSource,
  text: &'a str,
  run_out: Option<&'a Path>,
}

impl RawConfigVisitor for Prepare<'_> {
  type Output = Result<PreparedCommand, Report>;

  fn visit<R: RawConfig>(self) -> Self::Output {
    prepare::<R>(self.command, self.source, self.text, self.run_out)
  }
}

macro_rules! impl_raw_config {
  ($raw:ty, $args:ty, $variant:ident, $kind:ident, [$($required:expr),* $(,)?]) => {
    impl RawConfig for $raw {
      type Args = $args;

      fn wrap(args: Self::Args) -> CommandArgs {
        CommandArgs::$variant(Box::new(args))
      }

      fn into_config(self) -> CommandConfig {
        CommandConfig::$variant { config: Box::new(self) }
      }

      fn set_run_outputs(&mut self, out_dir: &Path) {
        self.output.output_all = Some(out_dir.to_path_buf());
        self.output_selection = run_output_selection(CommandKind::$kind, &self.output_selection, &[$($required),*]);
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

impl_raw_config!(
  TreetimeTimetreeArgsRaw,
  TreetimeTimetreeArgs,
  Timetree,
  Timetree,
  [
    TimetreeOutputSelection::Auspice,
    TimetreeOutputSelection::Tracelog,
    TimetreeOutputSelection::ClockCsv,
  ]
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
  TreetimeHomoplasyArgsRaw,
  TreetimeHomoplasyArgs,
  Homoplasy,
  Homoplasy,
  [
    HomoplasyOutputSelection::Auspice,
    HomoplasyOutputSelection::HomoplasyStats,
  ]
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

fn prepare<R: RawConfig>(
  command: AppCommand,
  source: &ConfigSource,
  text: &str,
  run_out: Option<&Path>,
) -> Result<PreparedCommand, Report> {
  let schema = command_schema::<R>();
  let merged = load_config_document::<R>(source, text, None)?;
  check_command_config(source, &merged, &schema)?;
  let mut raw: R = from_json_value(&merged)?;
  let mut defaults = R::default();
  if let Some(out_dir) = run_out {
    raw.set_run_outputs(out_dir);
    defaults.set_run_outputs(out_dir);
  }
  let config = settings_map(&raw)?;
  check_cli_rules(command, source, &config)?;
  let changed_settings = changed_settings(&config, &settings_map(&defaults)?, schema.as_value())?;
  let args = R::Args::try_from(raw)?;
  Ok(PreparedCommand {
    config,
    changed_settings,
    args: R::wrap(args),
  })
}

#[cfg(feature = "clap")]
fn check_cli_rules(command: AppCommand, source: &ConfigSource, config: &Map<String, Value>) -> Result<(), Report> {
  render_and_bail(source, "invalid configuration", cli_rule_diagnostics(command, config)?)
}

#[cfg(not(feature = "clap"))]
fn check_cli_rules(_command: AppCommand, _source: &ConfigSource, _config: &Map<String, Value>) -> Result<(), Report> {
  Ok(())
}

fn settings_map<R: Serialize>(raw: &R) -> Result<Map<String, Value>, Report> {
  match to_json_value(&raw)? {
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
  let Some((prefix, suffix)) = CDS_PLACEHOLDERS
    .iter()
    .find_map(|placeholder| name.split_once(placeholder))
  else {
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
