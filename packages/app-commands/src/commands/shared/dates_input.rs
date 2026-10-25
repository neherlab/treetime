use crate::commands::shared::input_warnings::warn_duplicate_names;
use crate::commands::shared::metadata::MetadataInput;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use treetime::progress::{LogSink, RunWarningKind};
use treetime::progress_warn;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::NodeMap;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_io::dates_csv::{DateConstraint, metadata_read_file};
use treetime_utils::datetime::options::DateParserOptions;

const UNMATCHED_NAMES_SHOWN: usize = 10;

pub(crate) fn read_input_dates(
  dates: &DatesInput,
  graph: &Graph,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<NodeMap<DateConstraint>, Report> {
  let DatesInput {
    metadata,
    date_column,
    date_format,
  } = dates;
  let path = &metadata.path;
  let rows = metadata_read_file(
    path,
    &metadata.delimiters,
    &metadata.id_columns,
    None,
    date_column.as_deref(),
  )
  .and_then(|table| table.dates(&date_parser_options(Some(date_format.as_str()))))
  .wrap_err("When reading dates")?;

  let pairing = pair_by_name(graph.get_nodes().map(|node| node.key()), names, rows);
  warn_duplicate_names(
    log,
    RunWarningKind::DuplicateMetadataNames,
    &format!("The metadata '{}' has more than one row named", path.display()),
    "TreeTime uses the first row of each name.",
    &pairing.duplicate_entry_names,
  );
  if !pairing.unmatched.is_empty() {
    let shown = pairing
      .unmatched
      .iter()
      .take(UNMATCHED_NAMES_SHOWN)
      .map(|(name, _)| name)
      .join(", ");
    let suffix = if pairing.unmatched.len() > UNMATCHED_NAMES_SHOWN {
      "..."
    } else {
      ""
    };
    progress_warn!(
      log,
      "Date constraints found for {} names not present in tree: {shown}{suffix}",
      pairing.unmatched.len()
    );
  }

  Ok(
    pairing
      .by_node
      .into_iter()
      .filter_map(|(key, date)| Some((key, date?)))
      .collect(),
  )
}

pub fn date_parser_options(date_format: Option<&str>) -> DateParserOptions {
  DateParserOptions {
    date_format: date_format.map(ToOwned::to_owned),
    ..DateParserOptions::default()
  }
}

#[derive(Debug, Clone)]
pub(crate) struct DatesInput {
  pub(crate) metadata: MetadataInput,
  pub(crate) date_column: Option<String>,
  pub(crate) date_format: String,
}
