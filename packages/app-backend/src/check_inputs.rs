use crate::command::AppCommandExt;
use crate::input_detection::{
  DateForm, DateReading, DetectedSettings, ReadDates, ReadInputs, ReadMetadata, detected_settings,
};
use crate::results::year_date::YearDate;
use app_commands::command::AppCommand;
use app_commands::commands::shared::alphabet::AlphabetNameCli;
use app_commands::commands::shared::dates_input::date_parser_options;
use app_commands::commands::shared::tree_input::TreeDialectArg;
use app_commands::json_value::SparseConfig;
use chrono::Datelike;
use deser::adapters::DisplayFromStr;
use deser::{Deserialize, Serialize};
use eyre::{Report, WrapErr};
use itertools::{Itertools, izip};
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::cmp::Reverse;
use std::collections::BTreeSet;
use std::ops::RangeInclusive;
use std::path::{Path, PathBuf};
use strum_macros::Display;
use treetime::alphabet::alphabet::Alphabet;
use treetime_io::csv::{DELIMITED_EXTENSIONS, default_metadata_delimiters, default_name_candidates};
use treetime_io::dates_csv::{DateConstraint, DateValue, MetadataTable, metadata_read_file};
use treetime_io::fasta::{FASTA_EXTENSIONS, fasta_read_file};
use treetime_io::nwk::NewickDialect;
use treetime_io::tree::{TREE_EXTENSIONS, tree_read_file};
use treetime_schema::skip_serializing_optionals;
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::{DATE_FORMATS, parse_date_with_format};
use treetime_utils::io::compression::COMPRESSION_EXTENSIONS;
use treetime_utils::io::json::from_json_value;

const ROUND_DAYS: [u32; 2] = [1, 15];

const NUMERIC_DATE_YEARS: RangeInclusive<f64> = 1000.0..=3000.0;

/// An input file of a command, as the form asks for it.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize)]
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CommandInput {
  pub kind: InputKind,
  pub need: InputNeed,
}

impl CommandInput {
  pub const fn new(kind: InputKind, need: InputNeed) -> Self {
    Self { kind, need }
  }
}

/// How much a run of the app needs an input file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize)]
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
  let command = request.command;
  let config = &request.config.0;
  let request = InputSettings::of(request.command, &request.config)?;
  let mut facts = InputFacts::default();
  let mut problem = |input: InputKind, alignment: Option<usize>, report: &Report| {
    facts.problems.push(InputProblem {
      input,
      alignment,
      message: format!("{report:#}"),
    });
  };

  let newick_dialect = request.newick_dialect.dialect();
  let tree = request
    .tree
    .as_deref()
    .map(|path| read_tree(path, newick_dialect))
    .transpose();
  let tree = tree.unwrap_or_else(|report| {
    problem(InputKind::Tree, None, &report);
    None
  });

  let alphabet_name = request.alphabet.unwrap_or_default();
  let alphabet = Alphabet::new(alphabet_name.into())?;
  let alignments = request
    .alignment
    .iter()
    .enumerate()
    .filter_map(|(index, path)| {
      let number = index + 1;
      match read_sequence_names(path, number, &alphabet) {
        Ok(read) => Some(read),
        Err(report) => {
          problem(InputKind::Alignment, Some(number), &report);
          None
        },
      }
    })
    .collect_vec();

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
  let date_options = date_parser_options(request.date_format.as_deref());
  let metadata = request
    .metadata
    .as_deref()
    .map(|path| {
      read_metadata(
        path,
        &delimiters,
        &id_candidates,
        request.date_column.as_deref(),
        &date_options,
      )
    })
    .transpose()
    .unwrap_or_else(|report| {
      problem(InputKind::Metadata, None, &report);
      None
    });

  let (metadata, metadata_names, date_problem, read_metadata) = match metadata {
    Some(read) => (Some(read.facts), Some(read.names), read.date_problem, Some(read.read)),
    None => (None, None, None, None),
  };
  if let Some(report) = date_problem {
    problem(InputKind::Metadata, None, &report);
  }

  let alignment_lengths = (!alignments.is_empty()
    && alignments.len() == request.alignment.len()
    && alignments
      .iter()
      .all(|(alignment, _)| alignment.min_length == alignment.max_length))
  .then(|| {
    alignments
      .iter()
      .map(|(alignment, _)| alignment.max_length)
      .collect_vec()
  });
  facts.detected = detected_settings(&ReadInputs {
    command,
    config,
    newick_dialect: tree.as_ref().map(|_| newick_dialect),
    metadata: read_metadata,
    alphabet: (!alignments.is_empty()).then_some(alphabet_name),
    alignment_lengths,
    sequence_length: request.sequence_length,
  })?;

  let tips = tree.as_ref().map(|(_, tips)| tips);
  facts.alignments = alignments
    .into_iter()
    .map(|(alignment, names)| AlignmentFacts {
      tips_without_sequence: tips.map(|tips| missing(tips, &names)),
      ..alignment
    })
    .collect();
  if let Some((tree, tips)) = tree {
    facts.tips_without_metadata = metadata_names.as_ref().map(|names| missing(&tips, names));
    facts.tree = Some(tree);
  }
  facts.metadata = metadata;
  Ok(facts)
}

/// A command configuration whose input files to inspect before a run.
#[derive(Clone, Debug, JsonSchema, Deserialize)]
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
  newick_dialect: TreeDialectArg,
  metadata: Option<PathBuf>,
  alignment: Vec<PathBuf>,
  metadata_id_columns: Vec<String>,
  metadata_delimiters: Vec<char>,
  date_column: Option<String>,
  date_format: Option<String>,
  alphabet: Option<AlphabetNameCli>,
  sequence_length: Option<usize>,
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
  /// Facts about each alignment file that could be read, in the order of the alignment setting.
  pub alignments: Vec<AlignmentFacts>,
  /// Facts about the metadata table, when it could be read.
  pub metadata: Option<MetadataFacts>,
  /// Tree tips without a metadata row; present when both the tree and the metadata could be read.
  pub tips_without_metadata: Option<Vec<String>>,
  /// What the readers used or found for each setting of how the input files are read, with the key of the setting
  /// that overrides it.
  pub detected: DetectedSettings,
  /// Inputs that could not be read, with the reason.
  pub problems: Vec<InputProblem>,
}

/// An input that could not be read.
#[derive(Clone, Debug, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct InputProblem {
  /// Input the problem concerns.
  pub input: InputKind,
  /// Number of the alignment file the problem concerns, for an alignment problem: its position in the alignment setting, from 1.
  pub alignment: Option<usize>,
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

fn read_tree(path: &Path, dialect: NewickDialect) -> Result<(TreeFacts, Vec<String>), Report> {
  let parsed = tree_read_file(path, dialect)?;
  let graph = &parsed.graph;
  let (unnamed_tips, named_tips): (Vec<_>, Vec<_>) = graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .partition(|key| parsed.unnamed_nodes.contains(key));
  let internal = graph.get_internal_nodes().collect_vec();
  let named_tips = named_tips
    .into_iter()
    .map(|key| parsed.names.pairing_name(key).to_owned())
    .collect_vec();
  let facts = TreeFacts {
    tips: graph.num_leaves(),
    internal_nodes: internal.len(),
    polytomies: internal.iter().filter(|node| node.outbound().len() > 2).count(),
    unnamed_tips: unnamed_tips.len(),
    duplicate_node_names: parsed.duplicate_names.clone(),
    tip_names: named_tips.iter().sorted().dedup().cloned().collect(),
  };
  Ok((facts, named_tips))
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
  /// Names of the named tips, sorted, each once: the names that settings such as `outgroup` can name.
  pub tip_names: Vec<String>,
}

fn read_sequence_names(
  path: &Path,
  number: usize,
  alphabet: &Alphabet,
) -> Result<(AlignmentFacts, BTreeSet<String>), Report> {
  let records = fasta_read_file(path, alphabet)?;
  let names = records.iter().map(|record| record.seq_name.clone()).collect_vec();
  let lengths = records
    .iter()
    .unique_by(|record| &record.seq_name)
    .map(|record| record.seq.len())
    .collect_vec();
  let facts = AlignmentFacts {
    number,
    path: path.to_path_buf(),
    sequences: records.len(),
    min_length: lengths.iter().copied().min().unwrap_or(0),
    max_length: lengths.iter().copied().max().unwrap_or(0),
    duplicate_names: duplicates(&names),
    tips_without_sequence: None,
  };
  Ok((facts, names.into_iter().collect()))
}

/// Facts about one alignment file.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct AlignmentFacts {
  /// Alignment number of the file: its position in the alignment setting, from 1.
  pub number: usize,
  /// Path of the file, as the alignment setting gives it.
  pub path: PathBuf,
  /// Number of sequences.
  pub sequences: usize,
  /// Length of the shortest sequence.
  pub min_length: usize,
  /// Length of the longest sequence; equal to `min_length` for an alignment.
  pub max_length: usize,
  /// Sequence names that occur more than once in the file.
  pub duplicate_names: Vec<String>,
  /// Tree tips without a sequence in the file; present when the tree could be read.
  pub tips_without_sequence: Option<Vec<String>>,
}

fn read_metadata(
  path: &Path,
  delimiters: &[char],
  id_candidates: &[String],
  date_column: Option<&str>,
  date_options: &DateParserOptions,
) -> Result<MetadataRead, Report> {
  metadata_read_file(path, delimiters, id_candidates, None, date_column)
    .map(|table| metadata_summary(table, date_options))
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    result_defaulted,
    reason = "a table without a date column has no date facts, which the summary shows as `None`"
  )
)]
pub(crate) fn metadata_summary(table: MetadataTable, date_options: &DateParserOptions) -> MetadataRead {
  let (dates, date_problem) = if table.date_column.is_ok() {
    match table.dates(date_options) {
      Ok(dates) => (Some(date_facts(&dates, date_options)), None),
      Err(report) => (None, Some(report)),
    }
  } else {
    (None, None)
  };
  let (dates, readings) = dates.unzip();
  let date_column = table.date_column.ok();
  MetadataRead {
    names: table.rows.iter().map(|row| row.name.clone()).collect(),
    read: ReadMetadata {
      delimiter: table.delimiter,
      id_column: table.id_column.clone(),
      date_column: date_column.clone(),
      dates: readings,
    },
    facts: MetadataFacts {
      rows: table.rows.len(),
      columns: table.columns,
      id_column: table.id_column,
      date_column,
      dates,
    },
    date_problem,
  }
}

pub(crate) struct MetadataRead {
  pub(crate) facts: MetadataFacts,
  pub(crate) names: BTreeSet<String>,
  pub(crate) read: ReadMetadata,
  pub(crate) date_problem: Option<Report>,
}

/// Facts about a metadata table.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
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

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    result_defaulted,
    reason = "a calendar date that the date parser read parses again in the format it matched"
  )
)]
fn date_facts(dates: &[(String, Option<DateConstraint>)], options: &DateParserOptions) -> (DateFacts, ReadDates) {
  let read = dates.iter().filter_map(|(_, date)| date.as_ref()).collect_vec();
  let forms = read.iter().map(|date| date_reading(date, options)).collect_vec();
  let exact_days = izip!(&read, &forms)
    .filter_map(|(date, (form, format))| {
      let format = format.as_deref().filter(|_| *form == DateForm::Calendar)?;
      parse_date_with_format(&date.raw, format, options)
        .ok()
        .map(|parsed| parsed.day())
    })
    .collect_vec();
  let readings = forms
    .into_iter()
    .counts()
    .into_iter()
    .map(|((form, format), count)| DateReading { form, format, count })
    .sorted_by_key(|reading| (Reverse(reading.count), reading.form, reading.format.clone()))
    .collect_vec();
  let facts = DateFacts {
    readable: read.len(),
    unreadable: dates
      .iter()
      .filter(|(_, date)| date.is_none())
      .map(|(name, _)| name.clone())
      .collect(),
    exact_days: exact_days.len(),
    on_day_1_or_15: exact_days.iter().filter(|day| ROUND_DAYS.contains(day)).count(),
    span: date_span(&read),
  };
  let summary = ReadDates {
    readings,
    read: read.len(),
    rows: dates.len(),
  };
  (facts, summary)
}

fn date_reading(date: &DateConstraint, options: &DateParserOptions) -> (DateForm, Option<String>) {
  match date.value {
    DateValue::Exact(_) => {
      let numeric = date.raw.parse::<f64>().ok().filter(|value| value.is_finite());
      match (numeric, calendar_format(&date.raw, options)) {
        (Some(value), _) if NUMERIC_DATE_YEARS.contains(&value) => (DateForm::Numeric, None),
        (_, Some(format)) => (DateForm::Calendar, Some(format)),
        (_, None) => (DateForm::Numeric, None),
      }
    },
    DateValue::Uncertain(_) => (DateForm::Uncertain, None),
    DateValue::Range(_) => (DateForm::Range, None),
  }
}

fn calendar_format(raw: &str, options: &DateParserOptions) -> Option<String> {
  options
    .date_format
    .iter()
    .map(String::as_str)
    .chain(DATE_FORMATS.iter().copied())
    .find(|format| parse_date_with_format(raw, format, options).is_ok())
    .map(str::to_owned)
}

fn date_span(dates: &[&DateConstraint]) -> Option<DateSpan> {
  let bounds = dates
    .iter()
    .map(|date| match date.value {
      DateValue::Exact(exact) => (exact.value, exact.value),
      DateValue::Uncertain(range) | DateValue::Range(range) => (range.start, range.end),
    })
    .collect_vec();
  let earliest = bounds.iter().map(|(earliest, _)| *earliest).reduce(f64::min)?;
  let latest = bounds.iter().map(|(_, latest)| *latest).reduce(f64::max)?;
  Some(DateSpan {
    earliest: YearDate::new(earliest),
    latest: YearDate::new(latest),
  })
}

/// Facts about the sampling dates of a metadata table.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct DateFacts {
  /// Samples with a date the date parser reads.
  pub readable: usize,
  /// Samples whose date cannot be read, by name.
  pub unreadable: Vec<String>,
  /// Samples whose date is a calendar day.
  pub exact_days: usize,
  /// Samples whose date is a calendar day on the 1st or 15th of a month, which often marks a date rounded to the month.
  pub on_day_1_or_15: usize,
  /// Earliest and latest sampling date that the parser read; a date range counts with its first and its last day.
  /// Absent when no date could be read.
  pub span: Option<DateSpan>,
}

/// The earliest and the latest of a set of dates.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize, Deserialize)]
pub struct DateSpan {
  /// The earliest date.
  pub earliest: YearDate,
  /// The latest date.
  pub latest: YearDate,
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
