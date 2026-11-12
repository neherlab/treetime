use app_commands::command::AppCommand;
use app_commands::commands::shared::alphabet::AlphabetNameCli;
use app_commands::commands::shared::partition::ConfigPartitions;
use app_commands::config::properties::setting_properties;
use deser::{Deserialize, Serialize};
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde_json::{Map, Value};
use std::collections::BTreeSet;
use treetime_io::nwk::NewickDialect;
use treetime_schema::skip_serializing_optionals;

const NEWICK_DIALECT: &str = "newick_dialect";
const METADATA_DELIMITERS: &str = "metadata_delimiters";
const METADATA_ID_COLUMNS: &str = "metadata_id_columns";
const DATE_COLUMN: &str = "date_column";
const DATE_FORMAT: &str = "date_format";
const ALPHABET: &str = "alphabet";
const SEQUENCE_LENGTH: &str = "sequence_length";
const ALIGNMENT: &str = "alignment";
const PARTITION: &str = "partition";

pub(crate) fn detected_settings(read: &ReadInputs<'_>) -> Result<DetectedSettings, Report> {
  let keys: BTreeSet<String> = setting_properties(read.command.config_schema().as_value())?
    .into_iter()
    .map(|property| property.key)
    .collect();
  let config = read.config;
  let offered = |key: &str| keys.contains(key);
  let metadata = read.metadata.as_ref().filter(|_| offered(METADATA_DELIMITERS));
  Ok(DetectedSettings {
    newick_dialect: read
      .newick_dialect
      .filter(|_| offered(NEWICK_DIALECT))
      .map(|dialect| DetectedDialect {
        dialect: dialect.to_string(),
        origin: SettingOrigin::chosen(config, NEWICK_DIALECT),
        setting: NEWICK_DIALECT.to_owned(),
      }),
    metadata_delimiter: metadata.map(|metadata| DetectedDelimiter {
      delimiter: metadata.delimiter.to_string(),
      origin: SettingOrigin::picked(config, METADATA_DELIMITERS),
      setting: METADATA_DELIMITERS.to_owned(),
    }),
    id_column: metadata.map(|metadata| DetectedColumn {
      column: metadata.id_column.clone(),
      origin: SettingOrigin::picked(config, METADATA_ID_COLUMNS),
      setting: METADATA_ID_COLUMNS.to_owned(),
    }),
    date_column: metadata
      .filter(|_| offered(DATE_COLUMN))
      .and_then(|metadata| metadata.date_column.clone())
      .map(|column| DetectedColumn {
        column,
        origin: SettingOrigin::found(config, DATE_COLUMN),
        setting: DATE_COLUMN.to_owned(),
      }),
    date_format: metadata
      .filter(|_| offered(DATE_FORMAT))
      .and_then(|metadata| metadata.dates.clone())
      .map(|dates| DetectedDateFormat {
        readings: dates.readings,
        read: dates.read,
        rows: dates.rows,
        origin: SettingOrigin::found(config, DATE_FORMAT),
        setting: DATE_FORMAT.to_owned(),
      }),
    alphabet: read
      .alphabet
      .filter(|_| offered(ALPHABET))
      .map(|alphabet| DetectedAlphabet {
        alphabet,
        origin: SettingOrigin::chosen(config, ALPHABET),
        setting: ALPHABET.to_owned(),
      }),
    sequence_length: sequence_length(read).filter(|fact| offered(&fact.setting)),
    partitions: partitions(read).filter(|_| offered(PARTITION)),
  })
}

/// What the input readers used or found for each setting of how the input files are read, each with the key of the
/// setting that overrides it. A fact is present when the command has its setting and the readers read the input it
/// concerns.
#[derive(Clone, Debug, Default, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct DetectedSettings {
  /// Newick dialect the tree was read in.
  pub newick_dialect: Option<DetectedDialect>,
  /// Delimiter of the columns of the metadata table.
  pub metadata_delimiter: Option<DetectedDelimiter>,
  /// Column of the metadata table that holds the sample names.
  pub id_column: Option<DetectedColumn>,
  /// Column of the metadata table that holds the sampling dates.
  pub date_column: Option<DetectedColumn>,
  /// Forms in which the date parser read the sampling dates, with the count of dates read.
  pub date_format: Option<DetectedDateFormat>,
  /// Alphabet the alignments were read in.
  pub alphabet: Option<DetectedAlphabet>,
  /// Number of sites of the sequences: the sites of every alignment file, or the `sequence_length` setting without an
  /// alignment.
  pub sequence_length: Option<DetectedLength>,
  /// Partitions of the sites of the alignment files.
  pub partitions: Option<DetectedPartitions>,
}

/// Where the value of a detected setting comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum SettingOrigin {
  /// The readers found the value in the input files, among the candidates of the setting.
  Detected,
  /// The configuration sets the value.
  Setting,
  /// The default of the setting, which the configuration does not change.
  Default,
}

impl SettingOrigin {
  fn chosen(config: &Map<String, Value>, key: &str) -> Self {
    if is_set(config, key) {
      Self::Setting
    } else {
      Self::Default
    }
  }

  fn found(config: &Map<String, Value>, key: &str) -> Self {
    if is_set(config, key) {
      Self::Setting
    } else {
      Self::Detected
    }
  }

  fn picked(config: &Map<String, Value>, key: &str) -> Self {
    match config.get(key) {
      Some(Value::Array(candidates)) if candidates.len() == 1 => Self::Setting,
      _ => Self::Detected,
    }
  }
}

/// The Newick dialect the tree was read in.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedDialect {
  /// The dialect as the setting spells it: `<structure>,<annotations>`, for example `classic,plain`.
  pub dialect: String,
  /// Where the dialect comes from.
  pub origin: SettingOrigin,
  /// Key of the setting that overrides the dialect.
  pub setting: String,
}

/// The delimiter of the columns of the metadata table.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedDelimiter {
  /// The delimiter character, for example a tab.
  pub delimiter: String,
  /// Where the delimiter comes from.
  pub origin: SettingOrigin,
  /// Key of the setting that overrides the delimiter.
  pub setting: String,
}

/// A column of the metadata table.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedColumn {
  /// Name of the column, as the table writes it.
  pub column: String,
  /// Where the column comes from.
  pub origin: SettingOrigin,
  /// Key of the setting that overrides the column.
  pub setting: String,
}

/// The forms in which the date parser read the sampling dates.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedDateFormat {
  /// Each form the parser read dates in, with its count, the most frequent form first.
  pub readings: Vec<DateReading>,
  /// Number of rows whose date the parser read.
  pub read: usize,
  /// Number of rows of the table.
  pub rows: usize,
  /// Where the date format comes from: `setting` when the configuration names a format that the parser tries first.
  pub origin: SettingOrigin,
  /// Key of the setting that overrides the date format.
  pub setting: String,
}

/// The dates that the parser read in one form.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct DateReading {
  /// Form of the dates.
  pub form: DateForm,
  /// The format of a calendar date, in the `strftime` notation of the `date_format` setting, for example `%Y-%m-%d`.
  pub format: Option<String>,
  /// Number of dates read in this form.
  pub count: usize,
  /// The form as a person reads it: the pattern of a calendar format, such as `YYYY-MM-DD` for `%Y-%m-%d`, or the
  /// name of the form, such as `decimal years`.
  pub text: String,
}

impl DateReading {
  pub fn new(form: DateForm, format: Option<String>, count: usize) -> Self {
    Self {
      text: reading_text(form, format.as_deref()),
      form,
      format,
      count,
    }
  }
}

/// Form of a sampling date.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, JsonSchema, Serialize, Deserialize)]
#[schemars(rename_all = "kebab-case")]
#[deser(rename_all = "kebab-case")]
pub enum DateForm {
  /// A numeric date, a decimal year such as `2015.3`.
  Numeric,
  /// A calendar day in the format that `format` names.
  Calendar,
  /// A date with unknown parts, such as `2015-XX-XX`, read as the range of days it allows.
  Uncertain,
  /// A range of two dates, such as `2015-01-01/2015-06-30`.
  Range,
}

/// The alphabet the alignments were read in.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedAlphabet {
  /// The alphabet.
  pub alphabet: AlphabetNameCli,
  /// Where the alphabet comes from.
  pub origin: SettingOrigin,
  /// Key of the setting that overrides the alphabet.
  pub setting: String,
}

/// The number of sites of the sequences.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedLength {
  /// Number of sites.
  pub length: usize,
  /// Where the number comes from: `detected` for the sites of the alignment files, `setting` for the
  /// `sequence_length` setting.
  pub origin: SettingOrigin,
  /// Key of the setting that changes the number: `alignment` for the sites of the alignment files.
  pub setting: String,
}

/// The partitions of the sites of the alignment files.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct DetectedPartitions {
  /// Each partition, in the order of the run.
  pub partitions: Vec<PartitionFacts>,
  /// Where the partitions come from: `default` for one partition per alignment file.
  pub origin: SettingOrigin,
  /// Key of the setting that overrides the partitions.
  pub setting: String,
}

/// One partition of the sites of an alignment file.
#[derive(Clone, Debug, PartialEq, Eq, JsonSchema, Serialize, Deserialize)]
pub struct PartitionFacts {
  /// Name of the partition in outputs and messages.
  pub label: String,
  /// Number of the alignment file of the partition, from 1.
  pub alignment: usize,
  /// Number of sites of the partition.
  pub sites: usize,
  /// Substitution model of the partition, as the `model` setting spells it.
  pub model: String,
}

pub(crate) struct ReadInputs<'a> {
  pub(crate) command: AppCommand,
  pub(crate) config: &'a Map<String, Value>,
  pub(crate) newick_dialect: Option<NewickDialect>,
  pub(crate) metadata: Option<ReadMetadata>,
  pub(crate) alphabet: Option<AlphabetNameCli>,
  pub(crate) alignment_lengths: Option<Vec<usize>>,
  pub(crate) sequence_length: Option<usize>,
}

pub(crate) struct ReadMetadata {
  pub(crate) delimiter: char,
  pub(crate) id_column: String,
  pub(crate) date_column: Option<String>,
  pub(crate) dates: Option<ReadDates>,
}

#[derive(Clone)]
pub(crate) struct ReadDates {
  pub(crate) readings: Vec<DateReading>,
  pub(crate) read: usize,
  pub(crate) rows: usize,
}

fn sequence_length(read: &ReadInputs<'_>) -> Option<DetectedLength> {
  match (&read.alignment_lengths, read.sequence_length) {
    (Some(lengths), _) => Some(DetectedLength {
      length: lengths.iter().sum(),
      origin: SettingOrigin::Detected,
      setting: ALIGNMENT.to_owned(),
    }),
    (None, Some(length)) => Some(DetectedLength {
      length,
      origin: SettingOrigin::Setting,
      setting: SEQUENCE_LENGTH.to_owned(),
    }),
    (None, None) => None,
  }
}

#[cfg_attr(
  dylint_lib = "custom",
  expect(
    error_dropped_by_pattern,
    reason = "the run checks report partitions that do not resolve, with the settings they concern"
  )
)]
fn partitions(read: &ReadInputs<'_>) -> Option<DetectedPartitions> {
  let lengths = read.alignment_lengths.as_ref()?;
  let Ok(plans) = ConfigPartitions::from_config(read.config).and_then(|partitions| partitions.resolve(lengths)) else {
    return None;
  };
  Some(DetectedPartitions {
    partitions: plans
      .iter()
      .map(|plan| PartitionFacts {
        label: plan.layout.label(),
        alignment: plan.layout.alignment,
        sites: plan.layout.columns.len(),
        model: plan.model.name().to_string().to_lowercase(),
      })
      .collect_vec(),
    origin: SettingOrigin::chosen(read.config, PARTITION),
    setting: PARTITION.to_owned(),
  })
}

fn is_set(config: &Map<String, Value>, key: &str) -> bool {
  config.get(key).is_some_and(|value| match value {
    Value::Null => false,
    Value::Array(items) => !items.is_empty(),
    _ => true,
  })
}

const FORMAT_PATTERNS: [(&str, &str); 10] = [
  ("%Y", "YYYY"),
  ("%y", "YY"),
  ("%m", "MM"),
  ("%d", "DD"),
  ("%e", "D"),
  ("%B", "MMMM"),
  ("%b", "MMM"),
  ("%h", "MMM"),
  ("%j", "DDD"),
  ("%%", "%"),
];

fn reading_text(form: DateForm, format: Option<&str>) -> String {
  match (form, format) {
    (DateForm::Calendar, Some(format)) => format_pattern(format),
    (DateForm::Calendar, None) => "calendar dates".to_owned(),
    (DateForm::Numeric, _) => "decimal years".to_owned(),
    (DateForm::Uncertain, _) => "uncertain dates".to_owned(),
    (DateForm::Range, _) => "date ranges".to_owned(),
  }
}

fn format_pattern(format: &str) -> String {
  let mut pattern = String::new();
  let mut rest = format;
  while !rest.is_empty() {
    if let Some((after, text)) = FORMAT_PATTERNS
      .iter()
      .find_map(|(specifier, text)| Some((rest.strip_prefix(specifier)?, text)))
    {
      pattern.push_str(text);
      rest = after;
    } else {
      let mut chars = rest.chars();
      pattern.extend(chars.next());
      rest = chars.as_str();
    }
  }
  pattern
}
