use chrono::Datelike;
use csv::{ReaderBuilder, Trim};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;
use std::path::{Path, PathBuf};
use strum_macros::Display;
use treetime::alphabet::alphabet::Alphabet;
use treetime_io::csv::{default_metadata_delimiters, default_name_candidates};
use treetime_io::dates_csv::{DateValue, DatesMap, read_dates};
use treetime_io::fasta::read_many_fasta_path;
use treetime_io::nwk::nwk_read_file;
use treetime_utils::datetime::options::DateParserOptions;
use treetime_utils::datetime::parse_date::parse_date;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::{make_error, make_internal_report};

const DATE_COLUMN_CANDIDATE: &str = "date";
const ROUND_DAYS: [u32; 2] = [1, 15];

/// Input files to inspect before a run, with the metadata settings of the command.
#[derive(Clone, Debug, Default, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct CheckInputsRequest {
  /// Newick tree.
  #[schemars(extend("x-path" = "input"))]
  pub tree: Option<PathBuf>,
  /// Metadata table with one row per sample.
  #[schemars(extend("x-path" = "input"))]
  pub metadata: Option<PathBuf>,
  /// FASTA alignment files.
  #[schemars(extend("x-path" = "input"))]
  pub alignment: Vec<PathBuf>,
  /// Candidate names of the metadata column that holds the sample names; the command defaults when empty.
  pub metadata_id_columns: Vec<String>,
  /// Candidate field delimiters of the metadata table; the command defaults when empty.
  pub metadata_delimiters: Vec<char>,
  /// Name of the metadata column that holds the sampling dates; detected when absent.
  pub date_column: Option<String>,
}

/// Facts about the input files of a run, read with the readers the commands use.
#[derive(Clone, Debug, Default, Serialize, Deserialize, JsonSchema)]
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct TreeFacts {
  /// Number of tips.
  pub tips: usize,
  /// Number of internal nodes, the root included.
  pub internal_nodes: usize,
  /// Number of internal nodes with more than two children.
  pub polytomies: usize,
  /// Number of tips without a name.
  pub unnamed_tips: usize,
  /// Tip names that occur more than once.
  pub duplicate_tip_names: Vec<String>,
}

/// Facts about an alignment.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
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
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct MetadataFacts {
  /// Number of data rows.
  pub rows: usize,
  /// Column names, in file order.
  pub columns: Vec<String>,
  /// Column that holds the sample names, when one was found.
  pub id_column: Option<String>,
  /// Column that holds the sampling dates, when one was found.
  pub date_column: Option<String>,
  /// Facts about the sampling dates, when the date column could be read.
  pub dates: Option<DateFacts>,
}

/// Facts about the sampling dates of a metadata table.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
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
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
pub struct InputProblem {
  /// Input the problem concerns.
  pub input: InputKind,
  /// The error, with its causes.
  pub message: String,
}

/// Kind of input file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema, Display)]
#[serde(rename_all = "kebab-case")]
#[strum(serialize_all = "kebab-case")]
pub enum InputKind {
  Tree,
  Metadata,
  Alignment,
}

pub fn check_inputs(request: &CheckInputsRequest) -> InputFacts {
  let mut facts = InputFacts::default();
  let mut problem = |input: InputKind, report: &Report| {
    facts.problems.push(InputProblem {
      input,
      message: format!("{report:#}"),
    });
  };

  let tree = request.tree.as_deref().map(read_tree).transpose();
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
    .map(|path| read_metadata(path, &delimiters, &id_candidates, request.date_column.as_ref()))
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
  facts
}

fn read_tree(path: &Path) -> Result<(TreeFacts, Vec<String>), Report> {
  let parsed = nwk_read_file(path)?;
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
    duplicate_tip_names: duplicates(&named_tips),
  };
  Ok((facts, named_tips))
}

fn read_sequence_names(paths: &[PathBuf]) -> Result<(AlignmentFacts, BTreeSet<String>), Report> {
  let records = read_many_fasta_path(paths, &Alphabet::default())?;
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

struct MetadataRead {
  facts: MetadataFacts,
  names: BTreeSet<String>,
  date_problem: Option<Report>,
}

fn read_metadata(
  path: &Path,
  delimiters: &[char],
  id_candidates: &[String],
  date_column: Option<&String>,
) -> Result<MetadataRead, Report> {
  let mut text = String::new();
  open_file_or_stdin(&Some(path))?
    .read_to_string(&mut text)
    .wrap_err_with(|| format!("When reading '{}'", path.display()))?;
  let table = parse_table(&text, path, delimiters, id_candidates)?;
  let id_column = find_column(&table.columns, id_candidates);
  let detected_date_column = match date_column {
    Some(column) => table.columns.iter().find(|name| *name == column).cloned(),
    None => find_column(&table.columns, &[DATE_COLUMN_CANDIDATE.to_owned()]),
  };
  let names = id_column
    .as_ref()
    .and_then(|column| table.columns.iter().position(|name| name == column))
    .map(|index| {
      table
        .rows
        .iter()
        .filter_map(|row| row.get(index).cloned())
        .collect::<BTreeSet<_>>()
    })
    .unwrap_or_default();

  let (dates, date_problem) = if detected_date_column.is_some() {
    match read_dates(path, delimiters, id_candidates, &None, &date_column.cloned()) {
      Ok(dates) => (Some(date_facts(&dates)), None),
      Err(report) => (None, Some(report)),
    }
  } else {
    (None, None)
  };

  Ok(MetadataRead {
    facts: MetadataFacts {
      rows: table.rows.len(),
      columns: table.columns,
      id_column,
      date_column: detected_date_column,
      dates,
    },
    names,
    date_problem,
  })
}

struct Table {
  columns: Vec<String>,
  rows: Vec<Vec<String>>,
}

fn parse_table(text: &str, path: &Path, delimiters: &[char], id_candidates: &[String]) -> Result<Table, Report> {
  let tables: Vec<Table> = delimiters
    .iter()
    .unique()
    .map(|delimiter| parse_with_delimiter(text, *delimiter))
    .try_collect()?;
  let with_ids = tables
    .iter()
    .positions(|table| find_column(&table.columns, id_candidates).is_some())
    .collect_vec();
  let chosen = match with_ids.as_slice() {
    [index] => Some(*index),
    _ => extension_delimiter(path)
      .and_then(|delimiter| delimiters.iter().unique().position(|candidate| *candidate == delimiter))
      .or_else(|| with_ids.first().copied()),
  };
  match chosen {
    Some(index) => tables
      .into_iter()
      .nth(index)
      .ok_or_else(|| make_internal_report!("metadata delimiter {index} has no parsed table")),
    None => make_error!(
      "no column of '{}' holds the sample names; looked for: {}",
      path.display(),
      id_candidates.join(", ")
    ),
  }
}

fn parse_with_delimiter(text: &str, delimiter: char) -> Result<Table, Report> {
  let delimiter = u8::try_from(u32::from(delimiter))
    .wrap_err_with(|| format!("Metadata delimiter {delimiter:?} must fit in one byte"))?;
  let mut reader = ReaderBuilder::new()
    .trim(Trim::All)
    .delimiter(delimiter)
    .flexible(true)
    .from_reader(text.as_bytes());
  let columns = reader
    .headers()?
    .iter()
    .map(|header| header.trim_start_matches('#').trim_end_matches('#').trim().to_owned())
    .collect_vec();
  let rows: Vec<Vec<String>> = reader
    .records()
    .map(|record| record.map(|record| record.iter().map(str::to_owned).collect_vec()))
    .try_collect()?;
  Ok(Table { columns, rows })
}

fn find_column(columns: &[String], candidates: &[String]) -> Option<String> {
  let candidates = candidates
    .iter()
    .map(|candidate| candidate.to_lowercase())
    .collect_vec();
  columns
    .iter()
    .find(|column| candidates.contains(&column.to_lowercase()))
    .cloned()
}

fn extension_delimiter(path: &Path) -> Option<char> {
  let name = path.file_name()?.to_string_lossy().to_lowercase();
  let name = [".gz", ".bz2", ".xz", ".zst"]
    .iter()
    .find_map(|suffix| name.strip_suffix(suffix))
    .unwrap_or(&name);
  match Path::new(name).extension()?.to_str()? {
    "csv" => Some(','),
    "tsv" => Some('\t'),
    "ssv" => Some(';'),
    _ => None,
  }
}

fn date_facts(dates: &DatesMap) -> DateFacts {
  let options = DateParserOptions::default();
  let exact_days: BTreeMap<&String, u32> = dates
    .iter()
    .filter_map(|(name, date)| {
      let date = date.as_ref()?;
      if !matches!(date.value, DateValue::Exact(_)) || date.raw.parse::<f64>().is_ok() {
        return None;
      }
      parse_date(&date.raw, &options).ok().map(|parsed| (name, parsed.day()))
    })
    .collect();
  DateFacts {
    readable: dates.values().filter(|date| date.is_some()).count(),
    unreadable: dates
      .iter()
      .filter(|(_, date)| date.is_none())
      .map(|(name, _)| name.clone())
      .collect(),
    exact_days: exact_days.len(),
    on_day_1_or_15: exact_days.values().filter(|day| ROUND_DAYS.contains(day)).count(),
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
