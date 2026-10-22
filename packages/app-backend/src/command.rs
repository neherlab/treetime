use crate::check_inputs::{CommandInput, InputKind, InputNeed};
use crate::run_checks::MissingSequences;
use crate::runs::errors::invalid;
use app_commands::command::{AppCommand, CommandArgs, RawConfig, settings_map};
use app_commands::command_config::CommandConfig;
use app_commands::commands::ancestral::args::TreetimeAncestralArgsRaw;
use app_commands::commands::clock::args::TreetimeClockArgsRaw;
use app_commands::commands::homoplasy::args::TreetimeHomoplasyArgsRaw;
use app_commands::commands::mugration::args::TreetimeMugrationArgsRaw;
use app_commands::commands::optimize::args::TreetimeOptimizeArgsRaw;
use app_commands::commands::prune::args::TreetimePruneArgsRaw;
use app_commands::commands::shared::output_args::{
  AncestralOutputSelection, ClockOutputSelection, HomoplasyOutputSelection, MugrationOutputSelection,
  OptimizeOutputSelection, PruneOutputSelection, TimetreeOutputSelection,
};
use app_commands::commands::timetree::args::TreetimeTimetreeArgsRaw;
use app_commands::config::cli_rules::cli_rule_diagnostics;
use app_commands::config::load::{
  check_command_config, invalid_setting_report, load_config_document, merge_value, read_settings,
};
use app_commands::config::properties::{PathRole, setting_properties};
use app_commands::config::schema::command_schema;
use app_commands::config::source::{ConfigSource, render_and_bail};
use app_commands::json_value::SparseConfig;
use app_output::output_plan::{CommandKind, OutputSelection, PathPlaceholder, ResolvedOutputs};
use deser::adapters::As;
use deser::{Deserialize, Serialize};
use deser_serde::Serde;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use strum::IntoEnumIterator;
use treetime::cancel::Cancel;
use treetime::progress::{LogSink, StageSink};
use treetime_schema::skip_serializing_optionals;
use treetime_utils::error::report_to_string;
use treetime_utils::io::json::{JsonPretty, json_write_str};
use treetime_utils::make_error;

const REQUEST_SOURCE: &str = "config.json";

impl AppCommandExt for AppCommand {
  fn inputs(self) -> &'static [CommandInput] {
    const TREE: CommandInput = CommandInput::new(InputKind::Tree, InputNeed::Required);
    const METADATA: CommandInput = CommandInput::new(InputKind::Metadata, InputNeed::Required);
    const ALIGNMENT: CommandInput = CommandInput::new(InputKind::Alignment, InputNeed::Required);
    const ALIGNMENT_RECOMMENDED: CommandInput = CommandInput::new(InputKind::Alignment, InputNeed::Recommended);
    const ALIGNMENT_OPTIONAL: CommandInput = CommandInput::new(InputKind::Alignment, InputNeed::Optional);
    match self {
      Self::Timetree => &[TREE, METADATA, ALIGNMENT_RECOMMENDED],
      Self::Clock | Self::Mugration => &[TREE, METADATA],
      Self::Ancestral | Self::Homoplasy | Self::Optimize => &[TREE, ALIGNMENT],
      Self::Prune => &[TREE, ALIGNMENT_OPTIONAL],
    }
  }

  fn uses_dates(self) -> bool {
    matches!(self, Self::Timetree | Self::Clock)
  }

  fn missing_sequences(self) -> MissingSequences {
    match self {
      Self::Ancestral | Self::Homoplasy => MissingSequences::LimitedUnlessIgnored,
      Self::Optimize | Self::Timetree | Self::Prune => MissingSequences::Limited,
      Self::Clock | Self::Mugration => MissingSequences::Unread,
    }
  }

  fn title(self) -> &'static str {
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

  fn description(self) -> &'static str {
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

  fn main_settings(self) -> &'static [&'static str] {
    match self {
      Self::Timetree => &[
        "clock_rate",
        "clock_std_dev",
        "confidence",
        "coalescent",
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

  fn config_over_defaults(self, settings: &Value) -> Result<Map<String, Value>, Report> {
    let mut config = Value::Object(self.default_config()?);
    if settings.is_object() {
      merge_value(&mut config, settings);
    }
    match config {
      Value::Object(config) => Ok(config),
      _ => make_error!("a command configuration must be a mapping of settings"),
    }
  }

  fn prepare_text(self, source_name: &str, text: &str) -> Result<PreparedCommand, Report> {
    prepare_source(self, source_name, text, None)
  }

  fn prepare_run(self, config: &Value, out_dir: &Path) -> Result<PreparedCommand, Report> {
    let mut config = config.clone();
    let Value::Object(settings) = &mut config else {
      return make_error!("a command configuration must be a mapping of settings");
    };
    self.remove_output_paths(settings)?;
    let text = json_write_str(&As::<_, Serde>::new(&config), JsonPretty(true))?;
    prepare_source(self, "config.json", &text, Some(out_dir))
  }

  fn remove_output_paths(self, settings: &mut Map<String, Value>) -> Result<(), Report> {
    for property in setting_properties(self.config_schema().as_value())? {
      if property.path_role == Some(PathRole::Output) {
        settings.remove(&property.key);
      }
    }
    Ok(())
  }
}

pub fn config_from_request(command: AppCommand, config: &SparseConfig) -> Result<CommandConfig, Report> {
  let value = Value::Object(config.0.clone());
  let source = ConfigSource::new(REQUEST_SOURCE, serde_json::to_string_pretty(&value)?);
  check_command_config(&source, &value, &command.config_schema())
    .map_err(|report| invalid(report_to_string(&report)))?;
  config_from_settings(command, &value).map_err(|report| invalid(report_to_string(&report)))
}

pub fn config_from_settings(command: AppCommand, settings: &Value) -> Result<CommandConfig, Report> {
  let settings = Value::Object(command.config_over_defaults(settings)?);
  command
    .command_config(&settings)
    .wrap_err_with(|| format!("When reading the configuration of a `{command}` run"))
}

pub fn run_command(
  args: &CommandArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<CommandOutcome, Report> {
  args.execute(cancel, stages, log)?;
  Ok(CommandOutcome {
    command: args.command(),
    output_files: written_files(&args.resolve_outputs()?, alignment_paths(args))?,
  })
}

pub trait AppCommandExt {
  fn inputs(self) -> &'static [CommandInput];

  fn uses_dates(self) -> bool;

  fn missing_sequences(self) -> MissingSequences;

  fn title(self) -> &'static str;

  fn description(self) -> &'static str;

  fn main_settings(self) -> &'static [&'static str];

  fn config_over_defaults(self, settings: &Value) -> Result<Map<String, Value>, Report>;

  fn prepare_text(self, source_name: &str, text: &str) -> Result<PreparedCommand, Report>;

  fn prepare_run(self, config: &Value, out_dir: &Path) -> Result<PreparedCommand, Report>;

  fn remove_output_paths(self, settings: &mut Map<String, Value>) -> Result<(), Report>;
}

fn prepare_source(
  command: AppCommand,
  source_name: &str,
  text: &str,
  run_out: Option<&Path>,
) -> Result<PreparedCommand, Report> {
  let source = ConfigSource::new(source_name, text);
  match command {
    AppCommand::Timetree => prepare::<TreetimeTimetreeArgsRaw>(command, &source, text, run_out),
    AppCommand::Optimize => prepare::<TreetimeOptimizeArgsRaw>(command, &source, text, run_out),
    AppCommand::Prune => prepare::<TreetimePruneArgsRaw>(command, &source, text, run_out),
    AppCommand::Ancestral => prepare::<TreetimeAncestralArgsRaw>(command, &source, text, run_out),
    AppCommand::Homoplasy => prepare::<TreetimeHomoplasyArgsRaw>(command, &source, text, run_out),
    AppCommand::Clock => prepare::<TreetimeClockArgsRaw>(command, &source, text, run_out),
    AppCommand::Mugration => prepare::<TreetimeMugrationArgsRaw>(command, &source, text, run_out),
  }
}

fn prepare<R: RunConfig>(
  command: AppCommand,
  source: &ConfigSource,
  text: &str,
  run_out: Option<&Path>,
) -> Result<PreparedCommand, Report> {
  let schema = command_schema::<R>();
  let merged = load_config_document::<R>(source, text, None)?;
  check_command_config(source, &merged, &schema)?;
  let mut raw: R = read_settings(&merged).map_err(|report| invalid_setting_report(source, report))?;
  let mut defaults = R::default();
  if let Some(out_dir) = run_out {
    raw.set_run_outputs(out_dir);
    defaults.set_run_outputs(out_dir);
  }
  let config = settings_map(&raw)?;
  render_and_bail(source, "invalid configuration", cli_rule_diagnostics(command, &config)?)?;
  let changed_settings = changed_settings(&config, &settings_map(&defaults)?, schema.as_value())?;
  let args = R::Args::try_from(raw)?;
  Ok(PreparedCommand {
    config,
    changed_settings,
    args: R::wrap(args),
  })
}

macro_rules! impl_run_config {
  ($raw:ty, $kind:ident, [$($required:expr),* $(,)?]) => {
    impl RunConfig for $raw {
      fn set_run_outputs(&mut self, out_dir: &Path) {
        self.output.output_all = Some(out_dir.to_path_buf());
        self.output_selection.output_selection =
          run_output_selection(CommandKind::$kind, &self.output_selection.output_selection, &[$($required),*]);
      }
    }
  };
}

impl_run_config!(
  TreetimeTimetreeArgsRaw,
  Timetree,
  [
    TimetreeOutputSelection::Auspice,
    TimetreeOutputSelection::Tracelog,
    TimetreeOutputSelection::ClockCsv,
    TimetreeOutputSelection::Gtr,
  ]
);
impl_run_config!(
  TreetimeOptimizeArgsRaw,
  Optimize,
  [OptimizeOutputSelection::Auspice, OptimizeOutputSelection::Gtr]
);
impl_run_config!(
  TreetimePruneArgsRaw,
  Prune,
  [PruneOutputSelection::Auspice, PruneOutputSelection::Gtr]
);
impl_run_config!(
  TreetimeAncestralArgsRaw,
  Ancestral,
  [AncestralOutputSelection::Auspice, AncestralOutputSelection::Gtr]
);
impl_run_config!(
  TreetimeHomoplasyArgsRaw,
  Homoplasy,
  [
    HomoplasyOutputSelection::Auspice,
    HomoplasyOutputSelection::HomoplasyStats,
  ]
);
impl_run_config!(TreetimeClockArgsRaw, Clock, [ClockOutputSelection::Auspice]);
impl_run_config!(
  TreetimeMugrationArgsRaw,
  Mugration,
  [MugrationOutputSelection::Auspice, MugrationOutputSelection::Gtr]
);

trait RunConfig: RawConfig {
  fn set_run_outputs(&mut self, out_dir: &Path);
}

fn changed_settings(
  config: &Map<String, Value>,
  defaults: &Map<String, Value>,
  schema: &Value,
) -> Result<Vec<String>, Report> {
  Ok(
    setting_properties(schema)?
      .into_iter()
      .filter(|property| property.path_role.is_none())
      .filter(|property| config.get(&property.key) != defaults.get(&property.key))
      .map(|property| property.key)
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

pub(crate) fn written_files(outputs: &ResolvedOutputs, alignments: &[PathBuf]) -> Result<Vec<OutputFile>, Report> {
  let partitions = partitions_by_path(outputs, alignments);
  let mut files = vec![];
  for (kind, paths) in outputs.paths_by_selection() {
    for planned in paths {
      let partition = partitions.get(&planned);
      files.extend(existing_files(&planned)?.into_iter().map(|path| OutputFile {
        path,
        kind,
        partition: partition.cloned(),
      }));
    }
  }
  Ok(files.into_iter().sorted().dedup_by(|a, b| a.path == b.path).collect())
}

fn partitions_by_path(outputs: &ResolvedOutputs, alignments: &[PathBuf]) -> BTreeMap<PathBuf, OutputPartition> {
  alignments
    .iter()
    .zip(outputs.partitions())
    .flat_map(|(alignment, number)| {
      let partition = OutputPartition {
        number,
        alignment: alignment.clone(),
      };
      let trees = outputs.partition_tree_outputs(number).into_values();
      let files = OutputSelection::iter()
        .filter_map(move |selection| outputs.partition_path(selection, number).map(Path::to_path_buf));
      trees.chain(files).map(move |path| (path, partition.clone()))
    })
    .collect()
}

fn alignment_paths(args: &CommandArgs) -> &[PathBuf] {
  match args {
    CommandArgs::Timetree(args) => &args.alignment.alignment,
    CommandArgs::Optimize(args) => &args.alignment.alignment,
    CommandArgs::Prune(args) => &args.alignment.alignment,
    CommandArgs::Ancestral(args) => &args.alignment.alignment,
    CommandArgs::Homoplasy(args) => &args.alignment.alignment,
    CommandArgs::Clock(_) | CommandArgs::Mugration(_) => &[],
  }
}

fn existing_files(planned: &Path) -> Result<Vec<PathBuf>, Report> {
  let name = planned
    .file_name()
    .map(|name| name.to_string_lossy().into_owned())
    .unwrap_or_default();
  let Some((prefix, suffix)) = PathPlaceholder::Cds
    .tokens()
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
    let path = entry
      .wrap_err_with(|| format!("When listing the outputs in '{}'", dir.display()))?
      .path();
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

pub struct PreparedCommand {
  pub config: Map<String, Value>,
  pub changed_settings: Vec<String>,
  pub args: CommandArgs,
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
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct OutputFile {
  /// Path of the file.
  pub path: PathBuf,
  /// Output selection that produced the file.
  pub kind: OutputSelection,
  /// Alignment file the output belongs to, for an output written once per alignment file.
  pub partition: Option<OutputPartition>,
}

/// The alignment file of an output written once per alignment file.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, JsonSchema, Serialize, Deserialize)]
pub struct OutputPartition {
  /// Partition number: the position of the file in the alignment setting, starting at 1.
  pub number: usize,
  /// Path of the alignment file.
  pub alignment: PathBuf,
}
