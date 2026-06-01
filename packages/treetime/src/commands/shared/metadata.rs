use serde::{Deserialize, Serialize};
use smart_default::SmartDefault;
use std::fmt::Debug;
use treetime_io::csv::default_name_candidates;

/// Default sampling-date string format (ISO 8601 calendar date), matching augur's `--date-format`.
pub const DEFAULT_DATE_FORMAT: &str = "%Y-%m-%d";

/// Metadata identity and delimiter options shared by every command that reads a metadata table
/// (`timetree`, `clock`, `mugration`).
///
/// `--metadata-id-columns` (alias `--name-column`) lists the candidate columns holding the taxon
/// identifier that links a metadata row to a tree tip; the first column present in the header wins.
/// Matching is case-insensitive (see `treetime-io` column detection). `--metadata-delimiters` lists
/// candidate field separators; the delimiter actually present in the file is used.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct MetadataIdArgs {
  /// Candidate column name(s) holding the taxon identifier that links metadata to tree tips
  ///
  /// The first listed column that is present in the header is used. Matching is case-insensitive.
  #[default(_code = "default_name_candidates()")]
  #[cfg_attr(
    feature = "clap",
    clap(
      long = "metadata-id-columns",
      visible_alias = "name-column",
      num_args = 1..,
      value_name = "COLUMN",
      default_values_t = default_name_candidates(),
    )
  )]
  pub metadata_id_columns: Vec<String>,

  /// Candidate field delimiter(s) for the metadata table
  ///
  /// The delimiter actually present in the file is used. Defaults to comma and tab.
  #[default(vec![',', '\t'])]
  #[cfg_attr(
    feature = "clap",
    clap(long = "metadata-delimiters", num_args = 1.., value_name = "CHAR", default_values_t = vec![',', '\t'])
  )]
  pub metadata_delimiters: Vec<char>,
}

/// Sampling-date column and parsing options shared by the time-aware commands (`timetree`, `clock`).
///
/// `--date-column` overrides auto-detection of the column holding sampling dates; when omitted, the
/// leftmost column whose name contains `date` (case-insensitive) is used. `--date-format` controls
/// parsing of string dates; numeric, ISO, and uncertain dates parse regardless.
#[derive(Debug, Clone, SmartDefault, Serialize, Deserialize)]
#[serde(default)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct DateColumnArgs {
  /// Label of the column to be used as sampling date (auto-detected when omitted)
  #[cfg_attr(feature = "clap", clap(long = "date-column", value_name = "COLUMN"))]
  pub date_column: Option<String>,

  /// Format used to parse string sampling dates (numeric, ISO, and uncertain dates parse regardless)
  #[default(DEFAULT_DATE_FORMAT.to_owned())]
  #[cfg_attr(feature = "clap", clap(long = "date-format", value_name = "FORMAT", default_value = DEFAULT_DATE_FORMAT))]
  pub date_format: String,
}
