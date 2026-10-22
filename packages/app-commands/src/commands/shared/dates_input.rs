use crate::commands::shared::input_warnings::warn_duplicate_names;
use crate::commands::shared::metadata::MetadataIdArgs;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::progress::{LogSink, RunWarningKind};
use treetime::progress_warn;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::pair_by_name;
use treetime_io::dates_csv::{DateConstraint, metadata_read_file};
use treetime_utils::datetime::options::DateParserOptions;

const UNMATCHED_NAMES_SHOWN: usize = 10;

pub fn date_parser_options(date_format: Option<&str>) -> DateParserOptions {
  DateParserOptions {
    date_format: date_format.map(ToOwned::to_owned),
    ..DateParserOptions::default()
  }
}

pub fn read_input_dates(
  path: &Path,
  metadata_id: &MetadataIdArgs,
  date_column: Option<&str>,
  date_format: &str,
  graph: &Graph,
  names: &NodeNames,
  log: &dyn LogSink,
) -> Result<BTreeMap<GraphNodeKey, DateConstraint>, Report> {
  let rows = metadata_read_file(
    path,
    &metadata_id.metadata_delimiters,
    &metadata_id.metadata_id_columns,
    None,
    date_column,
  )
  .and_then(|table| table.dates(&date_parser_options(Some(date_format))))
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
