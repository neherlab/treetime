use crate::command::AppCommand;
use crate::commands::shared::alignment::read_alignment;
use crate::commands::shared::tree_input::TreeDialectArg;
use crate::json_value::SparseConfig;
use chrono::Datelike;
use deser::adapters::DisplayFromStr;
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use strum_macros::Display;
use treetime::alphabet::alphabet::Alphabet;
use treetime_io::csv::{DELIMITED_EXTENSIONS, default_metadata_delimiters, default_name_candidates};
use treetime_io::dates_csv::{DateConstraint, DateValue, MetadataTable, metadata_read_file};
use treetime_io::fasta::FASTA_EXTENSIONS;
use treetime_io::nwk::NewickDialect;
use treetime_io::tree::{TREE_EXTENSIONS, tree_read_file};
use treetime_schema::skip_serializing_optionals;
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::parse_date;
use treetime_utils::io::compression::COMPRESSION_EXTENSIONS;
use treetime_utils::io::json::from_json_value;

const ROUND_DAYS: [u32; 2] = [1, 15];

/// A command configuration whose input files to inspect before a run.
#[derive(Clone, Debug, JsonSchema, Serialize, Deserialize)]
#[schemars(deny_unknown_fields)]
#[deser(deny_unknown_fields)]
pub struct CheckInputsRequest {
  /// Command the configuration is for.
  pub command: AppCommand,
  /// Settings of the command; the input files and the metadata settings are read from it, the rest is ignored.
  pub config: SparseConfig,
}

#[derive(Clone, Debug, Default, Deserialize)]
#[deser(default)]
struct InputSettings {
  tree: Option<PathBuf>,
  #[deser(as = DisplayFromStr)]
  tree_dialect: TreeDialectArg,
  metadata: Option<PathBuf>,
  alignment: Vec<PathBuf>,
  metadata_id_columns: Vec<String>,
  metadata_delimiters: Vec<char>,
  date_column: Option<String>,
}

impl InputSettings {
  fn of(command: AppCommand, config: &Map<String, Value>) -> Result<Self, Report> {
    let settings: Self = from_json_value(&Value::Object(config.clone())).wrap_err("When reading the input settings")?;
    let reads = |kind: InputKind| command.inputs().iter().any(|input| input.kind == kind);
    Ok(Self {
      tree: settings.tree.filter(|_| reads(InputKind::Tree)),
      metadata: settings.metadata.filter(|_| reads(InputKind::Metadata)),
      alignment: if reads(InputKind::Alignment) {
        settings.alignment
      } else {
        vec![]
      },
      ..settings
    })
  }
}

/// Facts about the input files of a run, read with the readers the commands use.
#[derive(Clone, Debug, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct InputFacts {
  /// Facts about the tree, when it could be read.
  pub tree: Option<TreeFacts>,
  /// Facts about the alignment, when it could be read.
  pub alignment: Option<AlignmentFacts>,
  /// Facts about the metadata table, when it could be read.
  pub metadata: Option<MetadataFacts>,
  /// Tree tips without a metadata row; present when both the tree and the metadata could be read.
  pub tips_without_metadata: Option<Vec<String>>,
  /// Tree tips without a sequence; present when both the tree and the alignment could be read.
  pub tips_without_sequence: Option<Vec<String>>,
  /// Inputs that could not be read, with the reason.
  pub problems: Vec<InputProblem>,
}

/// Facts about a tree.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct TreeFacts {
  /// Number of tips.
  pub tips: usize,
  /// Number of internal nodes, the root included.
  pub internal_nodes: usize,
  /// Number of internal nodes with more than two children.
  pub polytomies: usize,
  /// Number of tips without a name.
  pub unnamed_tips: usize,
  /// Node names that occur more than once, tips and internal nodes alike, sorted.
  pub duplicate_node_names: Vec<String>,
}

/// Facts about an alignment.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct AlignmentFacts {
  /// Number of sequences.
  pub sequences: usize,
  /// Length of the shortest sequence.
  pub min_length: usize,
  /// Length of the longest sequence; equal to `min_length` for an alignment.
  pub max_length: usize,
  /// Sequence names that occur more than once.
  pub duplicate_names: Vec<String>,
}

/// Facts about a metadata table.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct MetadataFacts {
  /// Number of data rows.
  pub rows: usize,
  /// Column names, in file order.
  pub columns: Vec<String>,
  /// Column that holds the sample names.
  pub id_column: String,
  /// Column that holds the sampling dates, when one was found.
  pub date_column: Option<String>,
  /// Facts about the sampling dates, when the date column could be read.
  pub dates: Option<DateFacts>,
}

/// Facts about the sampling dates of a metadata table.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DateFacts {
  /// Samples with a date the date parser reads.
  pub readable: usize,
  /// Samples whose date cannot be read, by name.
  pub unreadable: Vec<String>,
  /// Samples whose date is a calendar day.
  pub exact_days: usize,
  /// Samples whose date is a calendar day on the 1st or 15th of a month, which often marks a date rounded to the month.
  pub on_day_1_or_15: usize,
}

/// An input that could not be read.
#[derive(Clone, Debug, JsonSchema, Serialize, Deserialize)]
pub struct InputProblem {
  /// Input the problem concerns.
  pub input: InputKind,
  /// The error, with its causes.
  pub message: String,
}

/// Kind of input file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Display, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum InputKind {
  Tree,
  Metadata,
  Alignment,
}

impl InputKind {
  pub const fn setting(self) -> &'static str {
    match self {
      Self::Tree => "tree",
      Self::Metadata => "metadata",
      Self::Alignment => "alignment",
    }
  }

  pub const fn label(self) -> &'static str {
    match self {
      Self::Tree => "Tree",
      Self::Metadata => "Metadata",
      Self::Alignment => "Alignment",
    }
  }

  pub const fn formats(self) -> &'static str {
    match self {
      Self::Tree => "Newick or Nexus",
      Self::Metadata => "CSV, TSV or SSV table",
      Self::Alignment => "Aligned FASTA",
    }
  }

  pub fn extensions(self) -> Vec<String> {
    let formats: Vec<&str> = match self {
      Self::Tree => TREE_EXTENSIONS.to_vec(),
      Self::Metadata => DELIMITED_EXTENSIONS.iter().map(|(extension, _)| *extension).collect(),
      Self::Alignment => FASTA_EXTENSIONS.to_vec(),
    };
    formats
      .into_iter()
      .chain(COMPRESSION_EXTENSIONS)
      .map(str::to_owned)
      .collect()
  }
}

/// An input file an app command reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct CommandInput {
  /// Kind of the file; also the setting that names it.
  pub kind: InputKind,
  /// Whether a run of the app needs the file.
  pub need: InputNeed,
}

impl CommandInput {
  pub const fn new(kind: InputKind, need: InputNeed) -> Self {
    Self { kind, need }
  }
}

/// An input file of a command, as the form asks for it.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct InputSlot {
  /// Kind of the file; also the setting that names it.
  pub kind: InputKind,
  /// Whether a run of the app needs the file.
  pub need: InputNeed,
  /// Name of the input, for example `Tree`.
  pub label: String,
  /// File formats the readers accept, for example `Newick`.
  pub formats: String,
  /// File name extensions of the accepted formats, compressed forms included, without the dot.
  pub extensions: Vec<String>,
  /// Whether the setting takes a list of files.
  pub list: bool,
}

impl InputSlot {
  pub fn new(input: CommandInput, list: bool) -> Self {
    Self {
      kind: input.kind,
      need: input.need,
      label: input.kind.label().to_owned(),
      formats: input.kind.formats().to_owned(),
      extensions: input.kind.extensions(),
      list,
    }
  }
}

/// How much a run of the app needs an input file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum InputNeed {
  /// The run does not start without the file.
  Required,
  /// The run starts without the file, with less information.
  Recommended,
  /// The file adds information that the run can do without.
  Optional,
}

pub fn check_inputs(request: &CheckInputsRequest) -> Result<InputFacts, Report> {
  let request = InputSettings::of(request.command, &request.config)?;
  let mut facts = InputFacts::default();
  let mut problem = |input: InputKind, report: &Report| {
    facts.problems.push(InputProblem {
      input,
      message: format!("{report:#}"),
    });
  };

  let tree_dialect = request.tree_dialect.dialect();
  let tree = request
    .tree
    .as_deref()
    .map(|path| read_tree(path, tree_dialect))
    .transpose();
  let tree = tree.unwrap_or_else(|report| {
    problem(InputKind::Tree, &report);
    None
  });

  let sequences = (!request.alignment.is_empty())
    .then(|| read_sequence_names(&request.alignment))
    .transpose()
    .unwrap_or_else(|report| {
      problem(InputKind::Alignment, &report);
      None
    });

  let id_candidates = if request.metadata_id_columns.is_empty() {
    default_name_candidates()
  } else {
    request.metadata_id_columns.clone()
  };
  let delimiters = if request.metadata_delimiters.is_empty() {
    default_metadata_delimiters()
  } else {
    request.metadata_delimiters.clone()
  };
  let metadata = request
    .metadata
    .as_deref()
    .map(|path| read_metadata(path, &delimiters, &id_candidates, request.date_column.as_deref()))
    .transpose()
    .unwrap_or_else(|report| {
      problem(InputKind::Metadata, &report);
      None
    });

  let (metadata, metadata_names, date_problem) = match metadata {
    Some(read) => (Some(read.facts), Some(read.names), read.date_problem),
    None => (None, None, None),
  };
  if let Some(report) = date_problem {
    problem(InputKind::Metadata, &report);
  }

  if let Some((tree, tips)) = &tree {
    facts.tips_without_metadata = metadata_names.as_ref().map(|names| missing(tips, names));
    facts.tips_without_sequence = sequences.as_ref().map(|(_, names)| missing(tips, names));
    facts.tree = Some(tree.clone());
  }
  facts.alignment = sequences.map(|(alignment, _)| alignment);
  facts.metadata = metadata;
  Ok(facts)
}

fn read_tree(path: &Path, dialect: NewickDialect) -> Result<(TreeFacts, Vec<String>), Report> {
  let parsed = tree_read_file(path, dialect)?;
  let names = parsed.names();
  let graph = &parsed.graph;
  let tip_names = graph.get_leaves().map(|leaf| names[&leaf.key()].clone()).collect_vec();
  let internal = graph.get_internal_nodes().collect_vec();
  let named_tips = tip_names.iter().flatten().cloned().collect_vec();
  let facts = TreeFacts {
    tips: tip_names.len(),
    internal_nodes: internal.len(),
    polytomies: internal.iter().filter(|node| node.outbound().len() > 2).count(),
    unnamed_tips: tip_names.iter().filter(|name| name.is_none()).count(),
    duplicate_node_names: parsed.duplicate_names.clone(),
  };
  Ok((facts, named_tips))
}

fn read_sequence_names(paths: &[PathBuf]) -> Result<(AlignmentFacts, BTreeSet<String>), Report> {
  let records = read_alignment(paths, &Alphabet::default())?;
  let names = records.iter().map(|record| record.seq_name.clone()).collect_vec();
  let lengths = records.iter().map(|record| record.seq.len()).collect_vec();
  let facts = AlignmentFacts {
    sequences: records.len(),
    min_length: lengths.iter().copied().min().unwrap_or(0),
    max_length: lengths.iter().copied().max().unwrap_or(0),
    duplicate_names: duplicates(&names),
  };
  Ok((facts, names.into_iter().collect()))
}

pub(crate) struct MetadataRead {
  pub(crate) facts: MetadataFacts,
  pub(crate) names: BTreeSet<String>,
  pub(crate) date_problem: Option<Report>,
}

fn read_metadata(
  path: &Path,
  delimiters: &[char],
  id_candidates: &[String],
  date_column: Option<&str>,
) -> Result<MetadataRead, Report> {
  metadata_read_file(path, delimiters, id_candidates, None, date_column).map(metadata_summary)
}

pub(crate) fn metadata_summary(table: MetadataTable) -> MetadataRead {
  let (dates, date_problem) = if table.date_column.is_ok() {
    match table.dates() {
      Ok(dates) => (Some(date_facts(&dates)), None),
      Err(report) => (None, Some(report)),
    }
  } else {
    (None, None)
  };
  MetadataRead {
    names: table.rows.iter().map(|row| row.name.clone()).collect(),
    facts: MetadataFacts {
      rows: table.rows.len(),
      columns: table.columns,
      id_column: table.id_column,
      date_column: table.date_column.ok(),
      dates,
    },
    date_problem,
  }
}

fn date_facts(dates: &[(String, Option<DateConstraint>)]) -> DateFacts {
  let options = DateParserOptions::default();
  let exact_days = dates
    .iter()
    .filter_map(|(_, date)| {
      let date = date.as_ref()?;
      if !matches!(date.value, DateValue::Exact(_)) || date.raw.parse::<f64>().is_ok() {
        return None;
      }
      parse_date(&date.raw, &options).ok().map(|parsed| parsed.day())
    })
    .collect_vec();
  DateFacts {
    readable: dates.iter().filter(|(_, date)| date.is_some()).count(),
    unreadable: dates
      .iter()
      .filter(|(_, date)| date.is_none())
      .map(|(name, _)| name.clone())
      .collect(),
    exact_days: exact_days.len(),
    on_day_1_or_15: exact_days.iter().filter(|day| ROUND_DAYS.contains(day)).count(),
  }
}

fn missing(tips: &[String], names: &BTreeSet<String>) -> Vec<String> {
  tips.iter().filter(|tip| !names.contains(*tip)).cloned().collect()
}

fn duplicates(names: &[String]) -> Vec<String> {
  names
    .iter()
    .counts()
    .into_iter()
    .filter(|(_, count)| *count > 1)
    .map(|(name, _)| name.clone())
    .sorted()
    .collect()
}
